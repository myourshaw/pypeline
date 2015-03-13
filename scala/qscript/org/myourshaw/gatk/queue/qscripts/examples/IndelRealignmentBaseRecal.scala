package org.broadinstitute.sting.queue.qscripts.NL

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode

import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden
import java.io.FileNotFoundException
import scala.io.Source

import org.apache.commons.io.FilenameUtils

class IndelRealignmentBaseRecal extends QScript {
  qscript =>

  /****************************************************************************
  * Required Parameters
  *************************************************************************** */


  @Input(doc="input BAM file - or list of BAM files (can be specified more than once)", fullName="input", shortName="i", required=true)
  var input: List[File] = Nil

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: List[File] = Nil

  /****************************************************************************
  * Optional Parameters
  *************************************************************************** */

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=false)
  var indels: List[File] = Nil

  @Input(doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", fullName="project", shortName="p", required=false)
  var projectName: String = "project"

  @Input(doc="Output path for the processed BAM files.", fullName="output_directory", shortName="outputDir", required=false)
  var outputDir: String = ""

  @Input(doc="the -L interval string to be used by GATK - output bams at interval only", fullName="gatk_interval_string", shortName="L", required=false)
  var intervalString: String = ""

  @Input(doc="an intervals file to be used by GATK - output bams at intervals only", fullName="gatk_interval_file", shortName="intervals", required=false)
  var intervals: File = _

  @Input(doc="Cleaning model: KNOWNS_ONLY, USE_READS or USE_SW", fullName="clean_model", shortName="cm", required=false)
  var cleaningModel: String = "USE_READS"

  @Input(doc="Dont perform validation on the BAM files", fullName="no_validation", shortName="nv", required=false)
  var noValidation: Boolean = false

  @Input(fullName="maxReadsForRealignment", shortName="maxReads", doc="max reads allowed at an interval for realignment", required=false)
  var maxReads: Int = 20000   // 20000 is the default in IndelRealigner

  /****************************************************************************
   * Hidden Parameters
   *************************************************************************** */

  @Hidden
  @Input(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  @Hidden
  @Input(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""


  /****************************************************************************
  * Global Variables
  *************************************************************************** */

  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output

  var cleanModelEnum: ConsensusDeterminationModel = ConsensusDeterminationModel.USE_READS




  /****************************************************************************
  * Helper classes and methods
  *************************************************************************** */

  class ReadGroup (val id: String,
                   val lb: String,
                   val pl: String,
                   val pu: String,
                   val sm: String,
                   val cn: String,
                   val ds: String)
  {}

  def getSampleBamMap(bamFiles: List[File]): Map[String, File] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, File]
    for (bam <- bamFiles) {
      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if(!sampleTable.contains(sample))
          sampleTable(sample) = bam
        else {
          assert(sampleTable(sample) == bam, "Error!  " + sample + " is already associated with " 
                 + sampleTable(sample) + ", but I found the same sample in " + bam + ".")
        }
          
      }
    }
    return sampleTable.toMap
  }

  def getIndelCleaningModel(): ConsensusDeterminationModel = {
    if (cleaningModel == "KNOWNS_ONLY")
      ConsensusDeterminationModel.KNOWNS_ONLY
    else if (cleaningModel == "USE_SW")
      ConsensusDeterminationModel.USE_SW
    else
      ConsensusDeterminationModel.USE_READS
  }

  /****************************************************************************
  * Main script
  *************************************************************************** */


  def script = {
    // final output list of processed bam files
    var cohortList: List[File] = List()

    // sets the model for the Indel Realigner
    cleanModelEnum = getIndelCleaningModel()

    // keep a record of the number of contigs in the first bam file in the list
    val bams = 
      if (input.length == 1)
        QScriptUtils.createListFromFile(input(0))
      else
        input
                 
    if (nContigs < 0)
     nContigs = QScriptUtils.getNumberOfContigs(bams(0))

    // get sample bams
    val sampleBAMFiles: Map[String, File] = getSampleBamMap(bams)

    // if this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir, projectName + ".intervals")
    if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY)
      add(target(null, globalIntervals))

    // put each sample through the pipeline
    for ((sample, bam) <- sampleBAMFiles) {

      // BAM files generated by the pipeline
      val cleanedBam = swapExt(bam, ".bam", ".clean.bam")
      val dedupedBam = swapExt(bam, ".bam", ".clean.dedup.bam")
      val recalBam   = swapExt(bam, ".bam", ".clean.dedup.recal.bam")

      // Accessory files
      val targetIntervals = if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
      val metricsFile     = swapExt(bam, ".bam", ".metrics") // for dedup
      val preRecalFile    = swapExt(bam, ".bam", ".pre_recal.csv")
      val postRecalFile   = swapExt(bam, ".bam", ".post_recal.csv")
      val preOutPath      = swapExt(bam, ".bam", ".pre")
      val postOutPath     = swapExt(bam, ".bam", ".post")
      val preValidateLog  = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")


      // Validation is an optional step for the BAM file generated after
      // alignment and the final bam file of the pipeline.
      if (!noValidation) {
        add(validate(bam, preValidateLog),
            validate(recalBam, postValidateLog))
      }

      if (cleaningModel != ConsensusDeterminationModel.KNOWNS_ONLY)
        add(target(bam, targetIntervals))

      add(clean(bam, targetIntervals, cleanedBam),        // IndelRealignment
          dedup(cleanedBam, dedupedBam, metricsFile),
          cov(cleanedBam, preRecalFile),
          recal(cleanedBam, preRecalFile, recalBam),      // BaseQualRecalibration
          cov(recalBam, postRecalFile),
          analyzeCovariates(preRecalFile, preOutPath),
          analyzeCovariates(postRecalFile, postOutPath))

      cohortList :+= recalBam
    }

    // output a BAM list with all the processed per sample files
    val cohortFile = new File(qscript.outputDir, qscript.projectName + ".cohort.list")
    add(writeList(cohortList, cohortFile))
  }



  /****************************************************************************
  * Classes (GATK Walkers)
  *************************************************************************** */



  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 3
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
      this.maxRecordsInRam = 100000
  }

  case class target (inBam: File, outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (cleanModelEnum != ConsensusDeterminationModel.KNOWNS_ONLY)
      this.input_file = List(inBam)
    this.out = outIntervals
    this.mismatchFraction = 0.0
    if (!qscript.dbSNP.isEmpty)
      this.known ++= qscript.dbSNP
    if (indels != null)
      this.known ++= qscript.indels
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class clean (inBam: File, tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = List(inBam)
    this.targetIntervals = tIntervals
    this.out = outBam
    if (!qscript.dbSNP.isEmpty)
      this.known ++= qscript.dbSNP
    if (!qscript.indels.isEmpty)
      this.known ++= qscript.indels
    this.consensusDeterminationModel = cleanModelEnum
    this.compress = 1
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    if (!qscript.dbSNP.isEmpty)
      this.knownSites ++= qscript.dbSNP
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    if (!qscript.intervalString.isEmpty()) this.intervalsString ++= List(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    if (!qscript.intervalString.isEmpty()) this.intervalsString ++= List(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }



  /****************************************************************************
   * Classes (non-GATK programs)
   *************************************************************************** */


  case class analyzeCovariates (inRecalFile: File, outPath: File) extends AnalyzeCovariates {
    this.recal_file = inRecalFile
    this.output_dir = outPath.toString
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }

  case class dedup (inBam: File, outBam: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outBam
    this.metrics = metricsFile
    this.memoryLimit = 24
    this.analysisName = queueLogDir + outBam + ".dedup"
    this.jobName = queueLogDir + outBam + ".dedup"
  }

  case class joinBams (inBams: List[File], outBam: File) extends MergeSamFiles with ExternalCommonArgs {
    this.input = inBams
    this.output = outBam
    this.analysisName = queueLogDir + outBam + ".joinBams"
    this.jobName = queueLogDir + outBam + ".joinBams"
  }

  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam with ExternalCommonArgs {
    this.input :+= inSam
    this.output = outBam
    this.memoryLimit = 24
    this.sortOrder = sortOrderP
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class validate (inBam: File, outLog: File) extends ValidateSamFile with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outLog
    this.memoryLimit = 24
    this.REFERENCE_SEQUENCE = qscript.reference
    this.isIntermediate = false
    this.analysisName = queueLogDir + outLog + ".validate"
    this.jobName = queueLogDir + outLog + ".validate"
  }


  case class addReadGroup (inBam: File, outBam: File, readGroup: ReadGroup) extends AddOrReplaceReadGroups with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outBam
    this.RGID = readGroup.id
    this.RGCN = readGroup.cn
    this.RGDS = readGroup.ds
    this.RGLB = readGroup.lb
    this.RGPL = readGroup.pl
    this.RGPU = readGroup.pu
    this.RGSM = readGroup.sm
    this.analysisName = queueLogDir + outBam + ".rg"
    this.jobName = queueLogDir + outBam + ".rg"
  }

  case class revert (inBam: File, outBam: File, removeAlignmentInfo: Boolean) extends RevertSam with ExternalCommonArgs {
    this.output = outBam
    this.input :+= inBam
    this.removeAlignmentInformation = removeAlignmentInfo;
    this.sortOrder = if (removeAlignmentInfo) {SortOrder.queryname} else {SortOrder.coordinate}
    this.analysisName = queueLogDir + outBam + "revert"
    this.jobName = queueLogDir + outBam + ".revert"
  }

  case class convertToFastQ (inBam: File, outFQ: File) extends SamToFastq with ExternalCommonArgs {
    this.input :+= inBam
    this.fastq = outFQ
    this.analysisName = queueLogDir + outFQ + "convert_to_fastq"
    this.jobName = queueLogDir + outFQ + ".convert_to_fastq"
  }

  case class writeList(inBams: List[File], outBamList: File) extends ListWriterFunction {
    this.inputFiles = inBams
    this.listFile = outBamList
    this.analysisName = queueLogDir + outBamList + ".bamList"
    this.jobName = queueLogDir + outBamList + ".bamList"
  }
}

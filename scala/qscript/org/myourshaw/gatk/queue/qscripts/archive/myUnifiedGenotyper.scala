package org.myourshaw.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model

import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden

import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import org.apache.commons.io.FilenameUtils

class myUnifiedGenotyper extends QScript { qscript =>

  /****************************************************************************
   * Required Parameters.  All initialized to empty values.
   ****************************************************************************/

  @Input(fullName="input", shortName="i", doc="input BAM file - or list of BAM files", required=true)
  var input: File = _

  @Argument(fullName="vcfPrefix", shortName="V", doc="Prefix (path/name without extensions) of VCF files output by UnifiedGenotyper and subsequent analyses")
  var vcfPrefix: String = ""

  /****************************************************************************
   * Optional Parameters. Some initialized to defaults
   ****************************************************************************/

  @Input(fullName="reference", shortName="R", doc="Reference fasta file", required=false)
  var reference: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/1.2/b37/human_g1k_v37.fasta"

  @Input(fullName="dbsnp", shortName="D", doc="dbsnp ROD to use (must be in VCF format)", required=false)
  var dbSNP: File = "/share/apps/myourshaw/resources/dbsnp135/00-All.vcf"

  @Input(fullName="extra_indels", shortName="indels", doc="extra VCF files to use as reference indels for Indel Realignment", required=false)
  var indels: List[File] = Nil

  @Input(fullName="extra_realigned_bams", shortName="extra_bams", doc="Extra recalibrated bam files to use as additional data for Unified Genotyper and Variant Recalibration", required=false)
  var extra_bams: List[File] = Nil

  @Input(fullName="project", shortName="p", doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", required=false)
  var projectName: String = "project"

  @Input(fullName="output_directory", shortName="outputDir", doc="Output path for the processed BAM files.", required=false)
  var outputDir: String = ""

  @Input(fullName="gatk_interval_string", shortName="L", doc="the -L interval string to be used by GATK - output bams at interval only", required=false)
  var intervalString: String = ""

  @Input(fullName="gatk_interval_file", shortName="intervals", doc="an intervals file to be used by GATK - output bams at intervals only", required=false)
  var intervals: File = "/share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_donor6_acceptor27_cds_utr5_utr3.interval_list"

  @Input(fullName="clean_model", shortName="cm", doc="Cleaning model: KNOWNS_ONLY, USE_READS or USE_SW", required=false)
  var cleaningModel: String = "USE_READS"

  @Input(fullName="no_validation", shortName="nv", doc="Dont perform validation on the BAM files", required=false)
  var noValidation: Boolean = false

  @Argument(fullName="filterName", shortName="filters", doc="A optional list of filter names.", required=false)
  var filterName: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(fullName="filterExpression", shortName="filterExpression", doc="An optional list of filter expressions.", required=false)
  var filterExpression: List[String] = Nil

  /****************************************************************************
   * Hidden Parameters
   ****************************************************************************/
  @Hidden
  @Input(fullName="scatter_gather", shortName="sg", doc="How many ways to scatter/gather", required=false)
  var nContigs: Int = -1

  @Hidden
  @Input(doc="Define the default platform for Count Covariates -- useful for techdev purposes only.", fullName="default_platform", shortName="dp", required=false)
  var defaultPlatform: String = ""

  @Hidden
  @Input(doc="Run the pipeline in test mode only", fullName = "test_mode", shortName = "test", required=false)
  var testMode: Boolean = false


  /****************************************************************************
   * Global Variables
   ****************************************************************************/

  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output

  var cleanModelEnum: ConsensusDeterminationModel = ConsensusDeterminationModel.USE_READS

  /****************************************************************************
   * Helper classes and methods
   ****************************************************************************/

  class ReadGroup (val id: String,
                   val lb: String,
                   val pl: String,
                   val pu: String,
                   val sm: String,
                   val cn: String,
                   val ds: String)
  {}

  // Utility function to merge all bam files of similar samples. Generates one BAM file per sample.
  // It uses the sample information on the header of the input BAM files.
  //
  // Because the realignment only happens after these scripts are executed, in case you are using
  // bwa realignment, this function will operate over the original bam files and output over the
  // (to be realigned) bam files.
  def createSampleFiles(bamFiles: List[File], realignedBamFiles: List[File]): Map[String, List[File]] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, List[File]]
    val realignedIterator = realignedBamFiles.iterator
    for (bam <- bamFiles) {
      val rBam = realignedIterator.next  // advance to next element in the realignedBam list so they're in sync.

      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if (!sampleTable.contains(sample))
          sampleTable(sample) = List(rBam)
        else if ( !sampleTable(sample).contains(rBam))
          sampleTable(sample) :+= rBam
      }
    }
    return sampleTable.toMap
  }

  // Utility function to map samples to sample bam files.
  // It uses the sample information on the header of the input (recalibrated) sample BAM files.
  def mapSampleFiles(sampleFiles: List[File]): Map[String, List[File]] = {

    // Creating a table with SAMPLE information from each input BAM file
    val sampleTable = scala.collection.mutable.Map.empty[String, List[File]]
    for (bam <- sampleFiles) {
      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // only allow one sample per file. Bam files with multiple samples would require pre-processing of the file
      // with PrintReads to separate the samples. Tell user to do it himself!
      assert(!QScriptUtils.hasMultipleSamples(readGroups), "The pipeline requires that only one sample is present in a BAM file. Please separate the samples in " + bam)

      // Fill out the sample table with the readgroups in this file
      for (rg <- readGroups) {
        val sample = rg.getSample
        if (!sampleTable.contains(sample))
          sampleTable(sample) = List(bam)
        else if (!sampleTable(sample).contains(bam))
          sampleTable(sample) :+= bam
      }
    }
    return sampleTable.toMap
  }

  // Rebuilds the Read Group string to give BWA
  def addReadGroups(inBam: File, outBam: File, samReader: SAMFileReader) {
    val readGroups = samReader.getFileHeader.getReadGroups
    var index: Int = readGroups.length
    for (rg <- readGroups) {
      val intermediateInBam: File = if (index == readGroups.length) { inBam } else { swapExt(outBam, ".bam", index+1 + "-rg.bam") }
      val intermediateOutBam: File = if (index > 1) {swapExt(outBam, ".bam", index + "-rg.bam") } else { outBam}
      val readGroup = new ReadGroup(rg.getReadGroupId, rg.getLibrary, rg.getPlatform, rg.getPlatformUnit, rg.getSample, rg.getSequencingCenter, rg.getDescription)
      add(addReadGroup(intermediateInBam, intermediateOutBam, readGroup))
      index = index - 1
    }
  }

  def getIndelCleaningModel(): ConsensusDeterminationModel = {
    if (cleaningModel == "KNOWNS_ONLY")
      ConsensusDeterminationModel.KNOWNS_ONLY
    else if (cleaningModel == "USE_SW")
      ConsensusDeterminationModel.USE_SW
    else
      ConsensusDeterminationModel.USE_READS
  }

  def revertBams(bams: List[File], removeAlignmentInformation: Boolean): List[File] = {
    var revertedBAMList: List[File] = List()
    for (bam <- bams)
      revertedBAMList :+= revertBAM(bam, removeAlignmentInformation)
    return revertedBAMList
  }

  def revertBAM(bam: File, removeAlignmentInformation: Boolean): File = {
    val revertedBAM = swapExt(bam, ".bam", ".reverted.bam")
    add(revert(bam, revertedBAM, removeAlignmentInformation))
    return revertedBAM
  }

  /**
   * Exchanges the extension on a file without reverting to current directory, unlike the function in QScript.scala.
   * @param file File to look for the extension.
   * @param oldExtension Old extension to strip off, if present.
   * @param newExtension New extension to append.
   * @return new File with the new extension in the current directory.
   */
 override def swapExt(file: File, oldExtension: String, newExtension: String) =
    new File(file.getPath.stripSuffix(oldExtension) + newExtension)

 /****************************************************************************
   * Main script
   ****************************************************************************/

//  trait UnifiedGenotyperArguments extends CommandLineGATK {
//    this.reference_sequence = qscript.reference
//    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
//    // Set the memory limit to 2 gigabytes on each command.
//    this.memoryLimit = 6
//  }
def script() {
   // final output list of processed bam files
   var cohortList: List[File] = Nil

   // sets the model for the Indel Realigner
   cleanModelEnum = getIndelCleaningModel()

   // keep a record of the number of contigs in the first bam file in the list
   val sampleBams = QScriptUtils.createListFromFile(qscript.input)
   if (nContigs < 0)
     nContigs = QScriptUtils.getNumberOfContigs(sampleBams(0))

   // generate a BAM file per sample joining all per lane files if necessary
   val sampleBAMFiles: Map[String, List[File]] = mapSampleFiles(sampleBams)

   // VCF files generated by the pipeline
   val rawVcf =  qscript.vcfPrefix + ".ug.raw.vcf"

   // Accessory files
   val ugMetricsFile = swapExt(rawVcf, ".vcf", ".metrics")

   //val genotyper = ug(sampleBams, rawVcf, ugMetricsFile)
   //add(genotyper)
   add(ug(sampleBams, rawVcf, ugMetricsFile))
    /*
    evalUnfiltered.eval :+= genotyper.out
    evalUnfiltered.out = swapExt(genotyper.out, "vcf", "eval")

    variantFilter.variant = genotyper.out
    variantFilter.out = swapExt(qscript.bamFile, "bam", "filtered.vcf")
    variantFilter.filterName = filterName
    variantFilter.filterExpression = filterExpression

    evalFiltered.eval :+= variantFilter.out
    evalFiltered.out = swapExt(variantFilter.out, "vcf", "eval")

    add(genotyper, evalUnfiltered)
    // Only add variant filtration to the pipeline if filters were passed in
    if (filterNames.size > 0)
      add(variantFilter, evalFiltered)
    */

  cohortList ++= sampleBams
  // output a BAM list with all the processed per sample files
  val cohortFile = new File(qscript.outputDir + qscript.projectName + ".cohort.list")
  add(writeList(cohortList, cohortFile))
}

  /****************************************************************************
   * Classes (GATK Walkers)
   ****************************************************************************/

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
    if (!qscript.intervalString.isEmpty()) this.intervalsString ++= List(qscript.intervalString)
    else if (qscript.intervals != null) this.intervals :+= qscript.intervals
    this.logging_level = "INFO"
    this.memoryLimit = 6
    this.phone_home = GATKRunReport.PhoneHomeOption.STANDARD
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

  case class target (inBams: List[File], outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (cleanModelEnum != ConsensusDeterminationModel.KNOWNS_ONLY)
      this.input_file = inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known ++= List(qscript.dbSNP)
    if (indels != null)
      this.known ++= qscript.indels
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class clean (inBams: List[File], tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known ++= List(qscript.dbSNP)
    if (qscript.indels != null)
      this.known ++= qscript.indels
    this.consensusDeterminationModel = cleanModelEnum
    this.compress = 0
    this.noPGTag = qscript.testMode;
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    this.knownSites ++= List(qscript.dbSNP)
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
    this.no_pg_tag = qscript.testMode
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }

  case class ug (inBams: List[File], outVcf: File, metricsFile: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    this.input_file = inBams
    this.out = outVcf
    this.metrics_file = metricsFile
    this.dbsnp = qscript.dbSNP
    this.output_mode = OUTPUT_MODE.EMIT_VARIANTS_ONLY
    this.genotype_likelihoods_model = Model.BOTH
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 10.0
    this.dcov = 2500
    this.group = List("Standard")
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.num_threads = 8
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".ug"
    this.jobName = queueLogDir + outVcf + ".ug"
  }

  /****************************************************************************
   * Classes (non-GATK programs)
   ****************************************************************************/

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
    this.memoryLimit = 16
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
    this.sortOrder = sortOrderP
    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
  }

  case class validate (inBam: File, outLog: File) extends ValidateSamFile with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outLog
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

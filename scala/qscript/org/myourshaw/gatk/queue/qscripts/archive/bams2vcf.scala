/*
* Copyright (c) 2012-2013 Michael Yourshaw All Rights Reserved
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

//package org.myourshaw.sting.queue.qscripts

/**
 * Created by IntelliJ IDEA.
 * User: myourshaw
 * Date: 2012-02-15
 * Time: 16:13
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode
import org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model
import org.broadinstitute.sting.utils.variantcontext._
import org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode
import org.broadinstitute.sting.gatk.walkers.coverage.DoCOutputType

import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.function.ListWriterFunction
import org.broadinstitute.sting.commandline.Hidden

import org.broadinstitute.sting.gatk.phonehome.GATKRunReport
import org.apache.commons.io.FilenameUtils

class Bams2Vcf extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the Bams2Vcf.
  // 'qscript' is now the same as 'Bams2Vcf.this'
  qscript =>

  /****************************************************************************
   * Required Parameters.  All initialized to empty values.
   ****************************************************************************/

  //Required input is one or both of input (-I) and/or extra_recalibrated_bams (-i)
  @Input(fullName="input", shortName="I", doc="input aligned, unrecalibrated BAM file(s) with readgroups merged by library, duplicates marked, and de-duped library bam files merged by sample", required=false)
  var input: List[File] = Nil

  @Input(fullName="extra_recalibrated_bams", shortName="i", doc="input recalibrated sample bam file(s) to use as additional (or only) input to Unified Genotyper and subsequent analyses", required=false)
  var extra_recalibrated_bams: List[File] = Nil

  @Argument(fullName="vcfPrefix", shortName="V", doc="prefix (path/name without extensions) of VCF files output by UnifiedGenotyper and subsequent analyses", required=true)
  var vcf_prefix: String = ""

  /****************************************************************************
   * Optional Parameters. Some initialized to defaults
   ****************************************************************************/

  @Input(fullName="reference", shortName="R", doc="reference fasta file (default: /share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/human_g1k_v37_decoy.fasta)", required=false)
  var reference: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/human_g1k_v37_decoy.fasta"

  @Input(fullName="dbsnp", shortName="D", doc="dbsnp VCF file (default: /share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/dbsnp_137.b37.vcf)", required=false)
  var dbsnp: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/dbsnp_137.b37.vcf"

  @Input(fullName="bam_interval_string", shortName="bam_interval_string", doc="the -L interval string to be used by TableRecalibration - output recalibrated bams at interval only", required=false)
  var bam_interval_string: String = ""

  @Input(fullName="bam_interval_file", shortName="bam_interval_file", doc="an intervals file to be used by quality score recalibration - output variants at intervals only", required=false)
  var bam_interval_file: File = _

  @Input(fullName="doc_interval_file", shortName="doc_interval_file", doc="an intervals file to be used by DepthOfCoverage - count reads at intervals only (default: /share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_protein_coding_cds.interval_list)", required=false)
  var doc_interval_file: File = "/share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_protein_coding_cds.interval_list"

  @Input(fullName="doc_gene_list", shortName="doc_gene_list", doc="calculate the depth of coverage statistics over this list of genes (default: /share/apps/myourshaw/resources/refgene.b37.sorted.txt)", required=false)
  var doc_gene_list: File = "/share/apps/myourshaw/resources/refgene.b37.sorted.txt"

  @Input(fullName="ug_interval_file", shortName="L", doc="an intervals file to be used by Unified Genotyper - output variants at intervals only (default: /share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_donor6_acceptor27_cds_utr5_utr3.interval_list)", required=false)
  var ug_interval_file: File = "/share/apps/myourshaw/resources/intervals/EnsemblIntervals/intervals_Ensembl_genes_donor6_acceptor27_cds_utr5_utr3.interval_list"

  @Input(fullName="hapmap", shortName="hapmap", doc="a resource file to be used by VQSR as a training set - These high quality sites are used both to train the Gaussian mixture model and then again when choosing a LOD threshold based on sensitivity to truth sites (default: /scratch1/tmp/myourshaw/resources/gatk_resource_bundle/1.2/b37/hapmap_3.3.b37.sites.vcf)", required=false)
  var hapmap: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/hapmap_3.3.b37.sites.vcf"

  @Input(fullName="omni", shortName="omni", doc="a resource file to be used by VQSR as a training set - These polymorphic sites from the Omni genotyping array are used when training the model (default: /scratch1/tmp/myourshaw/resources/gatk_resource_bundle/1.2/b37/1000G_omni2.5.b37.sites.vcf)", required=false)
  var omni: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/1000G_omni2.5.b37.sites.vcf"

  @Input(fullName="mills", shortName="mills", doc="a resource file to be used by VQSR as a training set when modeling indels (default: /scratch1/tmp/myourshaw/resources/gatk_resource_bundle/1.2/b37/Mills_Devine_2hit.indels.b37.sites.vcf)", required=false)
  var mills: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf"

  @Input(fullName="known_indels", shortName="indels", doc="The current best set of known indels to be used for local realignment (note that we don't use dbSNP for this anymore)", required=false)
  var known_indels: List[File] = List("/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf", "/share/apps/myourshaw/resources/gatk_resource_bundle/2.5/b37/1000G_phase1.indels.b37.vcf")

  @Input(fullName="project", shortName="p", doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", required=false)
  var projectName: String = "project"

  @Input(fullName="output_directory", shortName="outputDir", doc="output path for the processed BAM files.", required=false)
  var outputDir: String = ""

  @Argument(fullName="clean_model", shortName="cm", doc="cleaning model: KNOWNS_ONLY, USE_READS or USE_SW (default: USE_READS)", required=false)
  var cleaningModel: String = "USE_READS"

  @Argument(fullName="no_validation", shortName="nv", doc="dont perform validation on the BAM files (default: false)", required=false)
  var noValidation: Boolean = false

  @Argument(fullName="is_whole_genome", shortName="whole_genome", doc="data is from exome capture libraries (default: false)", required=false)
  var is_whole_genome: Boolean = false

  @Argument(fullName="dedup_cleaned_bams", shortName="dedup", doc="dedup bam files after IndelRealigner (default: false)", required=false)
  var dedup_cleaned_bams: Boolean = false

  /****************************************************************************
   * Hidden Parameters
   ****************************************************************************/
  @Hidden
  @Argument(fullName="scatter_gather", shortName="sg", doc="how many ways to scatter/gather", required=false)
  var nContigs: Int = -1

//  @Hidden
//  @Input(fullName="default_platform", shortName="dp", doc="define the default platform for Count Covariates -- useful for techdev purposes only.", required=false)
//  var defaultPlatform: String = ""

//  @Hidden
//  @Input(fullName = "test_mode", shortName = "test", doc="run the pipeline in test mode only", required=false)
//  var testMode: Boolean = false

  /****************************************************************************
   * Global Variables
   ****************************************************************************/

  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output

  var cleanModelEnum: ConsensusDeterminationModel = ConsensusDeterminationModel.USE_READS

  var sampleCount = 0

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

def script() {

  assert(!(input.isEmpty && extra_recalibrated_bams.isEmpty), "No input bam files in input or extra_recalibrated_bams parameters")

   //count number of samples
   sampleCount = mapSampleFiles(input ::: extra_recalibrated_bams).keys.size
   //allBamFiles = input ::: extra_recalibrated_bams
   //// map each sample ID to bam file(s)
   //val mapAllBamFiles: Map[String, List[File]] = mapSampleFiles(allBamFiles)
   //sampleCount = mapAllBamFiles.keys.size

   var recalibratedBams = extra_recalibrated_bams.toSet
  
  // recalibrate unrecalibrated bam files, if any
  if (!input.isEmpty) {
    val unrecalibratedBams = input.toSet

    // set the model for the Indel Realigner
    cleanModelEnum = getIndelCleaningModel()

    // the number of contigs in the first bam file in the set (used to set scatterCount)
    if (nContigs < 0)
      nContigs = QScriptUtils.getNumberOfContigs(unrecalibratedBams.head)

    // if this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir + projectName + ".intervals")
    if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY)
      add(target(null, globalIntervals))

    // map each sample ID to bam file(s)
    val sampleBAMFiles: Map[String, List[File]] = mapSampleFiles(unrecalibratedBams.toList)

    for ((sample, bamList) <- sampleBAMFiles) {

      // BAM files generated by the pipeline
      val bam        = if (bamList.length == 1) bamList(0) else new File(qscript.projectName + "." + sample + ".bam")
      val cleanedBam = swapExt(bam, ".bam", ".realigned.bam")
      val dedupedBam = swapExt(bam, ".bam", ".realigned.dedup.bam")
      val recalBam   = if (dedup_cleaned_bams) swapExt(bam, ".bam", ".realigned.dedup.recalibrated.bam") else swapExt(bam, ".bam", ".realigned.recalibrated.bam")

      // Accessory files
      val targetIntervals = if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
      val dedupMetricsFile     = swapExt(bam, ".bam", ".dedup_metrics")
      val preRecalFile    = swapExt(bam, ".bam", ".pre_recal.csv")
      val postRecalFile   = swapExt(bam, ".bam", ".post_recal.csv")
      val preOutPath      = swapExt(bam, ".bam", ".pre")
      val postOutPath     = swapExt(bam, ".bam", ".post")
      val preValidateLog  = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")

      // Validation is an optional step for the BAM file generated after
      // alignment and the final bam file of the pipeline.
      if (!noValidation) {
        for (sampleFile <- bamList)
          add(validate(sampleFile, preValidateLog),
            validate(recalBam, postValidateLog))
      }

      if (cleaningModel != ConsensusDeterminationModel.KNOWNS_ONLY)
        add(target(bamList, targetIntervals))
      add(clean(bamList, targetIntervals, cleanedBam))

      if (dedup_cleaned_bams){
        add(dedup(cleanedBam, dedupedBam, dedupMetricsFile),
        cov(dedupedBam, preRecalFile),
        recal(dedupedBam, preRecalFile, recalBam))
      }
      else{
        add(cov(cleanedBam, preRecalFile),
        recal(cleanedBam, preRecalFile, recalBam))
      }

      add(cov(recalBam, postRecalFile),
      analyzeCovariates(preRecalFile, preOutPath),
      analyzeCovariates(postRecalFile, postOutPath))

      recalibratedBams += recalBam
    }
  }

  // union of new and extra recalibrated bams
  val recalibratedSampleBams = recalibratedBams.toList.sortWith(_.toString < _.toString)

  // variant genotyping and annotation on recalibrated input bams + previously recalibrated extra bams
  if (!recalibratedSampleBams.isEmpty){

    // the number of contigs in the first bam file in the set (used to set scatterCount)
    if (nContigs < 0)
      nContigs = QScriptUtils.getNumberOfContigs(recalibratedSampleBams.head)

    // VCF files generated by the pipeline
    val analysis_ready_vcf = qscript.vcf_prefix + ".analysis_ready.vcf"
    val annotated_vcf = qscript.vcf_prefix + ".annotated.vcf"

    // Accessory files
    val raw_vcf =  qscript.vcf_prefix + ".ug.raw.vcf"
    val doc_output_prefix =  qscript.vcf_prefix + ".coverage"
    val ug_metrics = swapExt(raw_vcf, ".vcf", ".metrics")
    val raw_snvs_vcf = vcf_prefix + ".ug.raw.snvs.vcf"
    val vqsr_recal = raw_snvs_vcf + ".recal"
    val tranches = raw_snvs_vcf + ".tranches"
    val vqsr_plots = raw_snvs_vcf + "plots.R"
    val recalibrated_snvs_vcf = qscript.vcf_prefix + ".recalibrated.snvs.vcf"
    val raw_indels_vcf = qscript.vcf_prefix  + ".ug.raw.indels.vcf"
    val filtered_indels_vcf = qscript.vcf_prefix + ".filtered.indels.vcf"
    val vqsr_indel_recal = raw_snvs_vcf + ".indel.recal"
    val indel_tranches = raw_snvs_vcf + ".indel.tranches"
    val vqsr_indel_plots = raw_snvs_vcf + ".indelplots.R"
    val recalibrated_indel_vcf = qscript.vcf_prefix + ".recalibrated.indels.vcf"
    val raw_vcf_eval = swapExt(raw_vcf, ".vcf",".VariantEval")
    val analysis_ready_vcf_eval = swapExt(analysis_ready_vcf, ".vcf",".VariantEval")

    val indel_filter_expression = if (sampleCount >= 10) List("QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0") else List("QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0")

    // Add command line functions to be run
    //scatter-gather may only work with depth by locus
    //add(depth(recalibratedSampleBams, doc_output_prefix))
    add(ug(recalibratedSampleBams, raw_vcf, ug_metrics))
    add(select(raw_vcf, raw_snvs_vcf, List(VariantContext.Type.SNP))
      ,select(raw_vcf, raw_indels_vcf, List(VariantContext.Type.INDEL)))
    add(vqsr(raw_snvs_vcf, vqsr_recal, tranches, vqsr_plots, Mode.SNP)
      ,applyvqsr(raw_snvs_vcf, vqsr_recal, tranches, recalibrated_snvs_vcf, Mode.SNP))
    if (!is_whole_genome){
      add(filter(raw_indels_vcf, filtered_indels_vcf, indel_filter_expression, List("GATKStandard")))
      add(combine(List(recalibrated_snvs_vcf, filtered_indels_vcf), analysis_ready_vcf))
    }
    else{
      add(vqsr(raw_indels_vcf, vqsr_indel_recal, indel_tranches, vqsr_indel_plots, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL)
        ,applyvqsr(raw_indels_vcf, vqsr_indel_recal, indel_tranches, recalibrated_indel_vcf, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL))
      add(combine(List(recalibrated_snvs_vcf, filtered_indels_vcf), analysis_ready_vcf))
    }
    add(eval(raw_vcf, raw_vcf_eval)
      ,eval(analysis_ready_vcf, analysis_ready_vcf_eval))
    add(annotate(analysis_ready_vcf, annotated_vcf))
    //TODO:
/*
    add(vax())
    add(pileup())
    add(zila_plots())
*/
  }

  // output a BAM list with all the processed per sample files
  val recalibratedSampleBamsListFile = new File(qscript.outputDir + qscript.projectName + ".recalibratedSampleBams.list")
  add(writeList(recalibratedSampleBams, recalibratedSampleBamsListFile))
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
    this.known = qscript.known_indels
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class clean (inBams: List[File], tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known = qscript.known_indels
    this.consensusDeterminationModel = cleanModelEnum
//    this.compress = 0
//    this.noPGTag = qscript.testMode;
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    this.knownSites ++= List(qscript.dbsnp)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
//    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    if (!qscript.bam_interval_string.isEmpty()) this.intervalsString ++= List(qscript.bam_interval_string)
    else if (qscript.bam_interval_file != null) this.intervals :+= qscript.bam_interval_file
    this.scatterCount = nContigs
    this.analysisName = queueLogDir + outRecalFile + ".covariates"
    this.jobName = queueLogDir + outRecalFile + ".covariates"
  }

  case class recal (inBam: File, inRecalFile: File, outBam: File) extends TableRecalibration with CommandLineGATKArgs {
    this.input_file :+= inBam
    this.recal_file = inRecalFile
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = outBam
    if (!qscript.bam_interval_string.isEmpty()) this.intervalsString ++= List(qscript.bam_interval_string)
    else if (qscript.bam_interval_file != null) this.intervals :+= qscript.bam_interval_file
//    this.no_pg_tag = qscript.testMode
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }
  
  case class depth (inBams: List[File], outPrefix: File) extends DepthOfCoverage with CommandLineGATKArgs {
    this.input_file = inBams
    this.out = outPrefix
    this.intervals = List(qscript.doc_interval_file)
    this.calculateCoverageOverGenes = qscript.doc_gene_list
    this.omitDepthOutputAtEachBase = true
    this.omitIntervalStatistics = true
    this.omitLocusTable = true
    this.dcov = 1000
    this.partitionType = List(DoCOutputType.Partition.sample, DoCOutputType.Partition.library, DoCOutputType.Partition.readgroup)
    this.summaryCoverageThreshold = List("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","40","50","60","70","80","90","100","110","120","130","140","150","160","170","180","190","200","300","400","500","600","700","800","900","1000")
    this.num_threads = 8
    //this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outPrefix + ".depth"
    this.jobName = queueLogDir + outPrefix + ".depth"
  }

  case class ug (inBams: List[File], outVcf: File, metricsFile: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    this.input_file = inBams
    this.out = outVcf
    this.metrics_file = metricsFile
    this.intervals = List(qscript.ug_interval_file)
    this.dbsnp = qscript.dbsnp
    this.read_filter = List("BadCigar")
    this.output_mode = OUTPUT_MODE.EMIT_VARIANTS_ONLY
    this.genotype_likelihoods_model = Model.BOTH
    this.standard_min_confidence_threshold_for_calling = 30.0
    this.standard_min_confidence_threshold_for_emitting = 10.0
    this.downsample_to_coverage = 250
    this.group = List("Standard")
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.num_threads = 8
    this.scatterCount = nContigs
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".ug"
    this.jobName = queueLogDir + outVcf + ".ug"
  }

  case class eval(inVcf: File, outEval: File) extends VariantEval with CommandLineGATKArgs {
    this.eval = List(inVcf)
    this.out = outEval
    this.dbsnp = qscript.dbsnp
    this.goldStandard = qscript.mills
    //Argument ST and ET has a bad value: The selected stratification Sample and evaluator VariantSummary are incompatible due to combinatorial memory requirements. Please disable one
    this.stratificationModule = List("Sample")
    this.doNotUseAllStandardModules = true
    this.evalModule = List("CompOverlap","CountVariants","IndelLengthHistogram","IndelSummary","MultiallelicSummary","TiTvVariantEvaluator","ValidationReport")
    this.isIntermediate = false
    this.analysisName = queueLogDir + outEval + ".eval"
    this.jobName = queueLogDir + outEval + ".eval"
  }

  case class select(inVcf: File, outVcf: File, seltype: List[VariantContext.Type]) extends SelectVariants with CommandLineGATKArgs {
    this.variant = inVcf
    this.out = outVcf
    this.selectTypeToInclude = seltype
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".select"
    this.jobName = queueLogDir + outVcf + ".select"
  }

  case class vqsr(inVcf: File, outRecal: File, outTranches: File, outPlots: File, recalmode: Mode) extends VariantRecalibrator with CommandLineGATKArgs {
    this.input :+= inVcf
    this.recal_file = outRecal
    this.tranches_file = outTranches
    this.rscript_file = outPlots
    this.mode = recalmode
    this.use_annotation ++= List("QD", "HaplotypeScore", "ReadPosRankSum", "FS")
    if (mode == Mode.SNP) {
      // http://www.broadinstitute.org/gsa/wiki/index.php/Managing_user_input#Example_usage_in_Queue_scripts
      this.resource :+= new TaggedFile(qscript.hapmap,"hapmap,VCF,known=false,training=true,truth=true,prior=15.0")
      this.resource :+= new TaggedFile(qscript.omni,"omni,VCF,known=false,training=true,truth=false,prior=12.0")
      this.resource :+= new TaggedFile(qscript.dbsnp,"dbsnp,VCF,known=true,training=false,truth=false,prior=6.0")
      this.use_annotation ++= List("MQRankSum", "MQ")
      if (is_whole_genome)
        this.use_annotation ++= List("DP")
    }
    else if (mode == Mode.INDEL)
      this.resource :+= new TaggedFile(qscript.mills,"mills,VCF,known=true,training=true,truth=true,prior=12.0")
    if (sampleCount >= 10)
      this.use_annotation ++= List("InbreedingCoeff")
    if (!is_whole_genome) {
      this.maxGaussians = if(sampleCount >= 30) 6 else 4
      if (sampleCount < 30)
        this.percentBadVariants = 0.05
    }
    this.isIntermediate = false
    this.analysisName = queueLogDir + outRecal + ".vqsr"
    this.jobName = queueLogDir + outRecal + ".vqsr"
  }

  case class applyvqsr(inVcf: File, inRecal: File, inTranches: File, outVcf: File, recalmode: Mode) extends ApplyRecalibration with CommandLineGATKArgs {
    this.input = List(inVcf)
    this.recal_file = inRecal
    this.tranches_file = inTranches
    this.out = outVcf
    this.ts_filter_level = 99.0
    this.mode = recalmode
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".applyvqsr"
    this.jobName = queueLogDir + outVcf + ".applyvqsr"
  }

  case class filter(inVcf: File, outVcf: File, filterNames: List[String], filterExpressions: List[String]) extends VariantFiltration with CommandLineGATKArgs {
    this.variant = inVcf
    this.out = outVcf
    this.filterName = filterNames
    this.filterExpression = filterExpressions
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".applyvqsr"
    this.jobName = queueLogDir + outVcf + ".applyvqsr"
  }

  case class combine(inVcfs: List[File], outVcf: File) extends CombineVariants with CommandLineGATKArgs {
    this.variant = inVcfs
    this.out = outVcf
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".combine"
    this.jobName = queueLogDir + outVcf + ".combine"
  }

  case class annotate(inVcf: File, outVcf: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.variant = inVcf
    this.out = outVcf
    this.dbsnp = qscript.dbsnp
    this.requireStrictAlleleMatch = true
    this.group = List("StandardAnnotation")
    //this.useAllAnnotations = true //requires snpEffFile
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".annotate"
    this.jobName = queueLogDir + outVcf + ".annotate"
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
  
  case class vax() {
    
  }

  case class pileup() {

  }

  case class zila_plots() {

  }
}

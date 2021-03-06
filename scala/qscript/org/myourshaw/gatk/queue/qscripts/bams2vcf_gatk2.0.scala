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
//import org.broadinstitute.sting.utils.variantcontext._
//import org.broadinstitute.sting.utils.variant._
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

  @Input(fullName="reference", shortName="R", doc="reference fasta file (default: /share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/human_g1k_v37_decoy.fasta)", required=false)
  var reference: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/human_g1k_v37_decoy.fasta"

  @Input(fullName="dbsnp", shortName="D", doc="dbsnp VCF file (default: /share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/dbsnp_138.b37.vcf)", required=false)
  var dbsnp: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/dbsnp_138.b37.vcf"

  //@Input(fullName="bam_interval_string", shortName="bam_interval_string", doc="the -L interval string to be used by TableRecalibration - output recalibrated bams at interval only", required=false)
  //var bam_interval_string: String = ""

  //@Input(fullName="bam_interval_file", shortName="bam_interval_file", doc="an intervals file to be used by quality score recalibration - output variants at intervals only", required=false)
  //var bam_interval_file: File = _

  @Input(fullName="doc_interval_file", shortName="doc_interval_file", doc="an intervals file to be used by DepthOfCoverage - count reads at intervals only (default: /share/apps/myourshaw/resources/intervals/merge/grand_unified_exome_ensembl_refgene_nextera_nimblegen_sureselect_b37.interval_list)", required=false)
  var doc_interval_file: File = "/share/apps/myourshaw/resources/intervals/merge/grand_unified_exome_ensembl_refgene_nextera_nimblegen_sureselect_b37.interval_list"

  @Input(fullName="doc_gene_list", shortName="doc_gene_list", doc="calculate the depth of coverage statistics over this list of genes (default: /share/apps/myourshaw/resources/refgene/refgene.b37.gatk_sorted_no-header.txt)", required=false)
  var doc_gene_list: File = "/share/apps/myourshaw/resources/refgene/refgene.b37.gatk_sorted_no-header.txt"

  //don't use something like the Nextera one, which chops off splice sites
  // @Input(fullName="ug_interval_file", shortName="L", doc="an intervals file to be used by Unified Genotyper - output variants at intervals only (default: /share/apps/myourshaw/resources/intervals/Nextera/NexteraRapidCapture_ExpandedExome_TargetedRegions.interval_list [same as TruSeq])", required=false)
  //var ug_interval_file: File = "/share/apps/myourshaw/resources/intervals/Nextera/NexteraRapidCapture_ExpandedExome_TargetedRegions.interval_list"
  @Input(fullName="ug_interval_file", shortName="L", doc="an intervals file to be used by Unified Genotyper - output variants at intervals only (default: /share/apps/myourshaw/resources/intervals/merge/grand_unified_exome_ensembl_refgene_nextera_nimblegen_sureselect_b37.interval_list)", required=false)
  var ug_interval_file: File = "/share/apps/myourshaw/resources/intervals/merge/grand_unified_exome_ensembl_refgene_nextera_nimblegen_sureselect_b37.interval_list"

  @Input(fullName="hapmap", shortName="hapmap", doc="a resource file to be used by VQSR as a training set - These high quality sites are used both to train the Gaussian mixture model and then again when choosing a LOD threshold based on sensitivity to truth sites (default: /scratch1/tmp/myourshaw/resources/gatk_resource_bundle/1.2/b37/hapmap_3.3.b37.vcf)", required=false)
  var hapmap: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/hapmap_3.3.b37.vcf"

  @Input(fullName="omni", shortName="omni", doc="a resource file to be used by VQSR as a training set - These polymorphic sites from the Omni genotyping array are used when training the model (default: /scratch1/tmp/myourshaw/resources/gatk_resource_bundle/1.2/b37/1000G_omni2.5.b37.vcf)", required=false)
  var omni: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/1000G_omni2.5.b37.vcf"

  @Input(fullName="onekg_snps", shortName="onekg_snps", doc="a resource file to be used by VQSR as a training set - the high confidence 1000 Genomes SNP calls (default: /share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/1000G_phase1.snps.high_confidence.b37.vcf)", required=false)
  var onekg_snps: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/1000G_phase1.snps.high_confidence.b37.vcf"

  @Input(fullName="onekg_indels", shortName="indels", doc="The 1000 Genomes indels (default: /share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/1000G_phase1.indels.b37.vcf)", required=false)
  var onekg_indels: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/1000G_phase1.indels.b37.vcf"

  @Input(fullName="mills", shortName="mills", doc="a resource file to be used by VQSR as a training set when modeling indels (default: /scratch1/tmp/myourshaw/resources/gatk_resource_bundle/1.2/b37/Mills_and_1000G_gold_standard.indels.b37.vcf)", required=false)
  var mills: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"

  @Argument(fullName="project", shortName="p", doc="the project name determines the final output (BAM file) base name. Example NA12878 yields NA12878.processed.bam", required=false)
  var projectName: String = "project"

  @Argument(fullName="output_directory", shortName="outputDir", doc="output path for the processed BAM files.", required=false)
  var outputDir: String = ""

  @Argument(fullName="clean_model", shortName="cm", doc="cleaning model: KNOWNS_ONLY, USE_READS or USE_SW (default: USE_READS)", required=false)
  var cleaningModel: String = "USE_READS"

  @Argument(fullName="no_validation", shortName="nv", doc="dont perform validation on the BAM files (default: false)", required=false)
  var noValidation: Boolean = false

  @Argument(fullName="is_whole_genome", shortName="whole_genome", doc="data is from exome capture libraries (default: false)", required=false)
  var is_whole_genome: Boolean = false

  @Argument(fullName="dedup_cleaned_bams", shortName="dedup", doc="dedup bam files after IndelRealigner (default: false)", required=false)
  var dedup_cleaned_bams: Boolean = false

  //ReduceReads is not in GATK 3.0
  @Argument(fullName="reduce_reads", shortName="reduce", doc="reduce bam files with ReduceReads (default: false)", required=false)
  var reduce_reads: Boolean = false

  @Argument(fullName="unified_genotyper", shortName="ug", doc="run Unified Genotyper (default: true)", required=false)
  var unified_genotyper: Boolean = true

  @Argument(fullName="haplotype_caller", shortName="hc", doc="run HaplotypeCaller (default: true)", required=false)
  var haplotype_caller: Boolean = true

  @Argument(fullName="depth_of_coverage", shortName="depth", doc="run DepthOfCoverage on bam files (default: false)", required=false)
  var depth_of_coverage: Boolean = false

  @Argument(fullName="available_nodes", shortName="nodes", doc="number of cluster nodes for scatter/gather (default: 93)", required=false)
  var available_nodes: Int = 93

  @Argument(fullName="scatter_count", shortName="scatter", doc="number of number of parallel jobs for scatter/gather (default: 0)", required=false)
  //if scatter_count argument is < 0, use -1 (no scatter)
  //if scatter_count argument is = 0 use max(nContigs, available_nodes)
  //if scatter_count argument is > 0 use scatter_count
  var scatter_count: Int = 0

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

  var cleanModelEnum: org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel = org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS

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

  def getIndelCleaningModel(): org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel = {
    if (cleaningModel == "KNOWNS_ONLY")
      org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.KNOWNS_ONLY
    else if (cleaningModel == "USE_SW")
      org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_SW
    else
      org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS
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
  val allBamFiles = input ::: extra_recalibrated_bams
  sampleCount = mapSampleFiles(allBamFiles).keys.size
  //allBamFiles = input ::: extra_recalibrated_bams
  //// map each sample ID to bam file(s)
  //val mapAllBamFiles: Map[String, List[File]] = mapSampleFiles(allBamFiles)
  //sampleCount = mapAllBamFiles.keys.size

  // the number of contigs in the first bam file in the set (may be used to set scatterCount)
  if (nContigs < 0) {
   nContigs = QScriptUtils.getNumberOfContigs(allBamFiles.head)
  }
  //if scatter_count argument is < 0, use -1 (no scatter)
  //if scatter_count argument is = 0 use max(nContigs, available_nodes)
  //if scatter_count argument is > 0 use scatter_count
  scatter_count = if (scatter_count > 0) scatter_count else if (scatter_count == 0 && nContigs > 1 && available_nodes > 1) nContigs.max(available_nodes) else -1

  var recalibratedBams = extra_recalibrated_bams.toSet

  // realign and recalibrate raw bam files, if any
  if (!input.isEmpty) {
    val unrecalibratedBams = input.toSet

    // set the model for the Indel Realigner
    cleanModelEnum = getIndelCleaningModel()

    // if this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir + qscript.projectName + ".intervals")
    if (cleaningModel == org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.KNOWNS_ONLY)
      add(realignerTargetCreator(null, globalIntervals))

    // map each sample ID to bam file(s)
    val sampleBAMFiles: Map[String, List[File]] = mapSampleFiles(unrecalibratedBams.toList)

    for ((sample, bamList) <- sampleBAMFiles) {

      // BAM files generated by the pipeline
      val bam = if (bamList.length == 1) bamList(0) else new File(qscript.projectName + "." + sample + ".bam")
      val realignedBam = swapExt(bam, ".bam", ".realigned.bam")
      val dedupedBam = swapExt(realignedBam, ".bam", ".dedup.bam")
      val recalibratedBam   = if (dedup_cleaned_bams) swapExt(bam, ".bam", ".realigned.dedup.recalibrated.bam") else swapExt(bam, ".bam", ".realigned.recalibrated.bam")
      //val reducedBam = swapExt(recalibratedBam, ".bam", ".reduced.bam")

      // Accessory files
      val targetIntervals = if (cleaningModel == org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.KNOWNS_ONLY) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
      val dedupMetricsFile = swapExt(bam, ".bam", ".dedup_metrics")
      val recal_dataFile = swapExt(bam, ".bam", ".recal_data.grp")
      //val preRecalFile = swapExt(bam, ".bam", ".pre_recal.csv")
      //val postRecalFile = swapExt(bam, ".bam", ".post_recal.csv")
      //val preOutPath = swapExt(bam, ".bam", ".pre")
      //val postOutPath = swapExt(bam, ".bam", ".post")
      val preValidateLog = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")

      // Validation is an optional step for the BAM file generated after
      // alignment and the final bam file of the pipeline.
      if (!noValidation) {
        for (sampleFile <- bamList)
          add(validate(sampleFile, preValidateLog),
            validate(recalibratedBam, postValidateLog))
      }

      //IndelRealigner
      if (cleaningModel != org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.KNOWNS_ONLY)
        add(realignerTargetCreator(bamList, targetIntervals))
      add(indelRealigner(bamList, targetIntervals, realignedBam))

      //BaseRecalibrator
      if (dedup_cleaned_bams){
        add(dedup(realignedBam, dedupedBam, dedupMetricsFile))
        //add(dedup(realignedBam, dedupedBam, dedupMetricsFile),
        //cov(dedupedBam, preRecalFile),
        //recal(dedupedBam, preRecalFile, recalibratedBam))
        add(baseRecalibrator(List(dedupedBam), recal_dataFile))
      }
      else{
        //add(cov(realignedBam, preRecalFile),
        //recal(realignedBam, preRecalFile, recalibratedBam))
        add(baseRecalibrator(List(realignedBam), recal_dataFile))
      }

      add(createRecalibrated(List(realignedBam), recal_dataFile, recalibratedBam))
      //add(cov(recalibratedBam, postRecalFile),
      //analyzeCovariates(preRecalFile, preOutPath),
      //analyzeCovariates(postRecalFile, postOutPath))

      recalibratedBams += recalibratedBam
    }
  }

  // union of new and extra recalibrated bams
  var recalibratedSampleBams: List[File] = List()

   //ReduceReads
   if (reduce_reads){
     var reducedBams: List[File] = List()
     for (recalBam <- recalibratedBams){
       val reducedBam = swapExt(recalBam, ".bam", ".reduced.bam")
       add(reduceReads(List(recalBam), reducedBam))
       reducedBams = reducedBam :: reducedBams
     }
     recalibratedSampleBams = reducedBams.sortWith(_.toString < _.toString)
   }
   else{
     recalibratedSampleBams = recalibratedBams.toList.sortWith(_.toString < _.toString)
   }

  // variant genotyping and annotation on recalibrated input bams + previously recalibrated extra bams
  if (!recalibratedSampleBams.isEmpty){

    //depth of coverage
    val doc_output_prefix =  qscript.vcf_prefix + ".coverage"

    //for UnifiedGenotyper
    //raw
    val raw_vcf_ug =  qscript.vcf_prefix + ".ug.raw.vcf"
    //eval raw
    val raw_vcf_eval_ug = swapExt(raw_vcf_ug, ".vcf",".VariantEval")
    //val ug_metrics = swapExt(raw_vcf_ug, ".vcf", ".metrics")
    //recalibrate snvs
    val raw_snvs_vcf_ug = vcf_prefix + ".ug.raw.snvs.vcf"
    val vqsr_snvs_recal_ug = raw_snvs_vcf_ug + ".recal"
    val vqsr_snvs_tranches_ug = raw_snvs_vcf_ug + ".tranches"
    val vqsr_snvs_plots_ug = raw_snvs_vcf_ug + "plots.R"
    val recalibrated_snvs_vcf_ug = qscript.vcf_prefix + ".ug.recalibrated.snvs.vcf"
    //recalibrate indels
    val raw_indels_vcf_ug = qscript.vcf_prefix  + ".ug.raw.indels.vcf"
    val vqsr_indels_recal_ug = raw_indels_vcf_ug + ".recal"
    val vqsr_indels_tranches_ug = raw_indels_vcf_ug + ".tranches"
    val vqsr_indels_plots_ug = raw_indels_vcf_ug + "plots.R"
    val recalibrated_indels_vcf_ug = qscript.vcf_prefix + ".ug.recalibrated.indels.vcf"
    //combine snvs & indels
    val analysis_ready_vcf_ug = qscript.vcf_prefix + ".ug.analysis_ready.vcf"
    //eval
    val analysis_ready_vcf_eval_ug = swapExt(analysis_ready_vcf_ug, ".vcf",".VariantEval")
    //annotate
    val annotated_vcf_ug = swapExt(analysis_ready_vcf_ug, ".vcf",".annotated.vcf")

    //for HaplotypeCaller
    //raw
    val raw_vcf_hc =  qscript.vcf_prefix + ".hc.raw.vcf"
    //eval raw
    val raw_vcf_eval_hc = swapExt(raw_vcf_hc, ".vcf",".VariantEval")
    //val hc_metrics = swapExt(raw_vcf_hc, ".vcf", ".metrics")
    //recalibrate snvs
    val raw_snvs_vcf_hc = vcf_prefix + ".hc.raw.snvs.vcf"
    val vqsr_snvs_recal_hc = raw_snvs_vcf_hc + ".recal"
    val vqsr_snvs_tranches_hc = raw_snvs_vcf_hc + ".tranches"
    val vqsr_snvs_plots_hc = raw_snvs_vcf_hc + "plots.R"
    val recalibrated_snvs_vcf_hc = qscript.vcf_prefix + ".hc.recalibrated.snvs.vcf"
    //recalibrate indels
    val raw_indels_vcf_hc = qscript.vcf_prefix  + ".hc.raw.indels.vcf"
    val vqsr_indels_recal_hc = raw_indels_vcf_hc + ".recal"
    val vqsr_indels_tranches_hc = raw_indels_vcf_hc + ".tranches"
    val vqsr_indels_plots_hc = raw_indels_vcf_hc + "plots.R"
    val recalibrated_indels_vcf_hc = qscript.vcf_prefix + ".hc.recalibrated.indels.vcf"
    //combine snvs & indels
    val analysis_ready_vcf_hc = qscript.vcf_prefix + ".hc.analysis_ready.vcf"
    //eval
    val analysis_ready_vcf_eval_hc = swapExt(analysis_ready_vcf_hc, ".vcf",".VariantEval")
    //annotate
    val annotated_vcf_hc = swapExt(analysis_ready_vcf_hc, ".vcf",".annotated.vcf")

    //val indel_filter_expression = if (sampleCount >= 10) List("QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0") else List("QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0")

    if(depth_of_coverage){
      //scatter-gather may only work with depth by locus
      add(depth(recalibratedSampleBams, doc_output_prefix))
    }

    //HaplotypeCaller
    if(haplotype_caller){
      add(haplotypeCaller(recalibratedSampleBams, raw_vcf_hc))
      //split snvs and indels into two files
      add(selectVariants(raw_vcf_hc, raw_snvs_vcf_hc, List(org.broadinstitute.variant.variantcontext.VariantContext.Type.SNP))
        ,selectVariants(raw_vcf_hc, raw_indels_vcf_hc, List(org.broadinstitute.variant.variantcontext.VariantContext.Type.INDEL)))
      //VQSR
      add(vqsr(raw_snvs_vcf_hc, vqsr_snvs_recal_hc, vqsr_snvs_tranches_hc, vqsr_snvs_plots_hc, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP, false)
        ,applyvqsr(raw_snvs_vcf_hc, vqsr_snvs_recal_hc, vqsr_snvs_tranches_hc, recalibrated_snvs_vcf_hc, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP))
      add(vqsr(raw_indels_vcf_hc, vqsr_indels_recal_hc, vqsr_indels_tranches_hc, vqsr_indels_plots_hc, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL, false)
        ,applyvqsr(raw_indels_vcf_hc, vqsr_indels_recal_hc, vqsr_indels_tranches_hc, recalibrated_indels_vcf_hc, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL))
      //combine snvs and indels
      add(combine(List(recalibrated_snvs_vcf_hc, recalibrated_indels_vcf_hc), analysis_ready_vcf_hc))
      //VariantEval
      add(variantEval(raw_vcf_hc, raw_vcf_eval_hc)
        ,variantEval(analysis_ready_vcf_hc, analysis_ready_vcf_eval_hc))
      //VariantAnnotator
      //as written this does nothing but zero out the FS annotation
      //maybe it will fix AD values messed up by CombineVariants?
      //add(annotate(analysis_ready_vcf_hc, annotated_vcf_hc))
    }

    //UnifiedGenotyper
    if (unified_genotyper){
      //add(unifiedGenotyper(recalibratedSampleBams, raw_vcf_ug, ug_metrics))
      add(unifiedGenotyper(recalibratedSampleBams, raw_vcf_ug))
      //split snvs and indels into two files
      add(selectVariants(raw_vcf_ug, raw_snvs_vcf_ug, List(org.broadinstitute.variant.variantcontext.VariantContext.Type.SNP))
        ,selectVariants(raw_vcf_ug, raw_indels_vcf_ug, List(org.broadinstitute.variant.variantcontext.VariantContext.Type.INDEL)))
      //VQSR
      add(vqsr(raw_snvs_vcf_ug, vqsr_snvs_recal_ug, vqsr_snvs_tranches_ug, vqsr_snvs_plots_ug, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP, true)
        ,applyvqsr(raw_snvs_vcf_ug, vqsr_snvs_recal_ug, vqsr_snvs_tranches_ug, recalibrated_snvs_vcf_ug, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP))
      add(vqsr(raw_indels_vcf_ug, vqsr_indels_recal_ug, vqsr_indels_tranches_ug, vqsr_indels_plots_ug, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL, true)
        ,applyvqsr(raw_indels_vcf_ug, vqsr_indels_recal_ug, vqsr_indels_tranches_ug, recalibrated_indels_vcf_ug, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL))
      /*
      if (!is_whole_genome){
        add(filter(raw_indels_vcf_ug, filtered_indels_vcf_ug, indel_filter_expression, List("GATKStandard")))
        add(combine(List(recalibrated_snvs_vcf_ug, filtered_indels_vcf_ug), analysis_ready_vcf_ug))
      }
      else{
        add(vqsr(raw_indels_vcf_ug, vqsr_indel_recal_ug, indel_tranches_ug, vqsr_indel_plots_ug, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL)
          ,applyvqsr(raw_indels_vcf_ug, vqsr_indel_recal_ug, indel_tranches_ug, recalibrated_indel_vcf_ug, org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL))
        add(combine(List(recalibrated_snvs_vcf_ug, filtered_indels_vcf_ug), analysis_ready_vcf_ug))
      }
      */
      //combine snvs and indels although AD values may me wrong
      add(combine(List(recalibrated_snvs_vcf_ug, recalibrated_indels_vcf_ug), analysis_ready_vcf_ug))
      //VariantEval
      add(variantEval(raw_vcf_ug, raw_vcf_eval_ug)
        ,variantEval(analysis_ready_vcf_ug, analysis_ready_vcf_eval_ug))
      //VariantAnnotator
      //as written this does nothing but zero out the FS annotation
      //maybe it will fix AD values messed up by CombineVariants?
      //add(annotate(analysis_ready_vcf_ug, annotated_vcf_ug))
    }

    //TODO:
/*
    add(vax())
    add(pileup())
*/
  }

  //output a list of bam files used in analysis
  val recalibrated_sample_bams_list = qscript.outputDir + qscript.projectName + ".recalibratedSampleBams.list"
  val recalibratedSampleBamsListFile = new File(recalibrated_sample_bams_list)
  add(writeList(recalibratedSampleBams, recalibratedSampleBamsListFile))
}

  /****************************************************************************
   * Classes (GATK Walkers)
   ****************************************************************************/

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 6
    this.isIntermediate = true
  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
    this.logging_level = "INFO"
    this.memoryLimit = 6
    this.phone_home = GATKRunReport.PhoneHomeOption.AWS
  }

  // General arguments to GATK walkers to be run on high memory node
  trait CommandLineGATKArgsHiMem extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
    this.logging_level = "INFO"
    this.memoryLimit = 20
    this.phone_home = GATKRunReport.PhoneHomeOption.AWS
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

//http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools

  case class realignerTargetCreator (inBams: List[File], outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (cleanModelEnum != org.broadinstitute.sting.gatk.walkers.indels.IndelRealigner.ConsensusDeterminationModel.KNOWNS_ONLY)
      this.input_file = inBams
    this.out = outIntervals
    this.mismatchFraction = 0.0
    this.known =List(qscript.mills, qscript.onekg_indels)
    this.scatterCount = qscript.scatter_count
    this.num_threads = 4 //maybe 6 or even 8?
    this.isIntermediate = false
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class indelRealigner (inBams: List[File], tIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known = List(qscript.mills, qscript.onekg_indels)
    this.consensusDeterminationModel = cleanModelEnum
//    this.compress = 0
//    this.noPGTag = qscript.testMode;
    this.scatterCount = qscript.scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class baseRecalibrator(inBams: List[File], recal_data: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.input_file = inBams
    this.knownSites ++= List(qscript.dbsnp, qscript.mills, qscript.onekg_indels)
    this.out = recal_data
    //this parameter seems to has been dropped
    //omit plot_pdf_file if using Queue scatter/gather
    //this.plot_pdf_file = recal_data+".plot.pdf"
    this.read_filter = List("BadCigar")
    //this should work but fails in gather combine tables
    //this.scatterCount = qscript.scatter_count
    //test this
    this.num_cpu_threads_per_data_thread = 8
    //this may be necessary to get plots, if plots are still possible
    //this.isIntermediate = true
    this.isIntermediate = false
    this.analysisName = queueLogDir + recal_data + ".baseRecalibrator"
    this.jobName = queueLogDir + recal_data + ".baseRecalibrator"
  }

  case class createRecalibrated(inBams: List[File], recal_data: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file = inBams
    this.BQSR = recal_data
    this.out = outBam
    //this.scatterCount = qscript.scatter_count
    this.num_cpu_threads_per_data_thread = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".createRecalibrated"
    this.jobName = queueLogDir + outBam + ".createRecalibrated"
  }

  case class reduceReads(inBams: List[File], outBam: File) extends ReduceReads with CommandLineGATKArgs {
    this.input_file = inBams
    this.out = outBam
    // with scatter/gather there is no guarantee that compressed read name uniqueness will be maintained
    this.dont_compress_read_names = true
    this.scatterCount = qscript.scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".reduceReads"
    this.jobName = queueLogDir + outBam + ".reduceReads"
  }

  /*
  case class cov (inBam: File, outRecalFile: File) extends CountCovariates with CommandLineGATKArgs {
    this.knownSites ++= List(qscript.dbsnp)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = outRecalFile
//    if (!defaultPlatform.isEmpty) this.default_platform = defaultPlatform
    if (!qscript.bam_interval_string.isEmpty()) this.intervalsString ++= List(qscript.bam_interval_string)
    else if (qscript.bam_interval_file != null) this.intervals :+= qscript.bam_interval_file
    this.scatterCount = qscript.scatter_count
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
    this.scatterCount = qscript.scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".recalibration"
    this.jobName = queueLogDir + outBam + ".recalibration"
  }
  */

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
    //this.scatterCount = qscript.scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outPrefix + ".depth"
    this.jobName = queueLogDir + outPrefix + ".depth"
  }

  //case class unifiedGenotyper (inBams: List[File], outVcf: File, metricsFile: File) extends UnifiedGenotyper with CommandLineGATKArgs {
  case class unifiedGenotyper (inBams: List[File], outVcf: File) extends UnifiedGenotyper with CommandLineGATKArgs {
    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    this.input_file = inBams
    this.out = outVcf
    //this.metrics_file = metricsFile
    this.intervals = List(qscript.ug_interval_file)
    this.dbsnp = qscript.dbsnp
    this.read_filter = List("BadCigar")
    this.output_mode = org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY
    this.genotype_likelihoods_model = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.BOTH
    this.standard_min_confidence_threshold_for_calling = 30.0
    this.standard_min_confidence_threshold_for_emitting = 10.0
    this.downsample_to_coverage = 250
    this.group = List("Standard")
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    //number of cores requested =  num_threads * num_cpu_threads_per_data_thread
    this.num_threads = 2
    //this.num_cpu_threads_per_data_thread = 12
    this.num_cpu_threads_per_data_thread = 4
    //each thread seems to consume ~2Gb
    this.jobNativeArgs  = Seq("-w n -l excl=true -l vf=2G")
    this.scatterCount = qscript.scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".ug"
    this.jobName = queueLogDir + outVcf + ".ug"
  }

  case class haplotypeCaller (inBams: List[File], outVcf: File) extends HaplotypeCaller with CommandLineGATKArgsHiMem {
    //the default memoryLimit of 6 may be OK for ReduceReads BAMs more may be required for uncompressed BAMs
    this.input_file = inBams
    this.out = outVcf
    this.intervals = List(qscript.ug_interval_file)
    this.dbsnp = qscript.dbsnp
    this.read_filter = List("BadCigar")
    this.standard_min_confidence_threshold_for_calling = 30.0
    this.standard_min_confidence_threshold_for_emitting = 10.0
    this.downsample_to_coverage = 250
    this.group = List("Standard")
    //##### ERROR MESSAGE: Invalid command line: Argument baq has a bad value: Walker cannot accept BAQ'd base qualities, and yet BAQ mode CALCULATE_AS_NECESSARY was requested.
    //this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    //##### ERROR MESSAGE: Invalid command line: Argument nt has a bad value: The analysis HaplotypeCaller currently does not support parallel execution with nt.  Please run your analysis without the nt option.
    //this.num_threads = 8
    this.jobNativeArgs  = Seq("-w n -l excl=true -l mem_free=20G")
    this.scatterCount = qscript.scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".hc"
    this.jobName = queueLogDir + outVcf + ".hc"}

  case class variantEval(inVcf: File, outEval: File) extends VariantEval with CommandLineGATKArgs {
    this.eval = List(inVcf)
    this.out = outEval
    this.dbsnp = qscript.dbsnp
    this.goldStandard = qscript.mills
    //Argument ST and ET has a bad value: The selected stratification Sample and evaluator VariantSummary are incompatible due to combinatorial memory requirements. Please disable one
    this.stratificationModule = List("Sample", "Filter")
    this.doNotUseAllStandardModules = true
    this.evalModule = List("CompOverlap","CountVariants","IndelLengthHistogram","IndelSummary","MultiallelicSummary","TiTvVariantEvaluator","ValidationReport")
    this.requireStrictAlleleMatch = true
    this.num_threads = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outEval + ".eval"
    this.jobName = queueLogDir + outEval + ".eval"
  }

  case class selectVariants(inVcf: File, outVcf: File, seltype: List[org.broadinstitute.variant.variantcontext.VariantContext.Type]) extends SelectVariants with CommandLineGATKArgs {
    this.variant = inVcf
    this.out = outVcf
    this.selectTypeToInclude = seltype
    this.num_threads = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".select"
    this.jobName = queueLogDir + outVcf + ".select"
  }

  case class vqsr(inVcf: File, outRecal: File, outTranches: File, outPlots: File, recalmode: org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode, input_from_ug: Boolean) extends VariantRecalibrator with CommandLineGATKArgs {
    this.input :+= inVcf
    this.recal_file = outRecal
    this.tranches_file = outTranches
    this.rscript_file = outPlots
    this.mode = recalmode
    this.percentBad = if(sampleCount >= 30) 0.01 else 0.05
    this.minNumBad = 1000
    if (is_whole_genome) {
      this.use_annotation ++= List("DP")
    }
    if (sampleCount >= 10){
      this.use_annotation ++= List("InbreedingCoeff")
    }
    if (mode == org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.SNP) {
      //For indels ##### ERROR MESSAGE: Bad input: Values for HaplotypeScore annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations. See http://gatkforums.broadinstitute.org/discussion/49/using-variant-annotator
      if (input_from_ug) {
        this.use_annotation ++= List("HaplotypeScore")
      }
      // http://www.broadinstitute.org/gsa/wiki/index.php/Managing_user_input#Example_usage_in_Queue_scripts
      this.resource :+= new TaggedFile(qscript.hapmap,"hapmap,VCF,known=false,training=true,truth=true,prior=15.0")
      this.resource :+= new TaggedFile(qscript.omni,"omni,VCF,known=false,training=true,truth=true,prior=12.0")
      this.resource :+= new TaggedFile(qscript.onekg_snps,"onekg_snps,VCF,known=false,training=true,truth=false,prior=10.0")
      this.resource :+= new TaggedFile(qscript.dbsnp,"dbsnp,VCF,known=true,training=false,truth=false,prior=2.0")
      this.use_annotation ++= List("QD", "MQRankSum", "MQ", "ReadPosRankSum", "FS")
      if (sampleCount < 30) {
        this.maxGaussians = 4
      }
    }
    else if (mode == org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode.INDEL) {
      this.resource :+= new TaggedFile(qscript.mills,"mills,VCF,known=false,training=true,truth=true,prior=12.0")
      this.resource :+= new TaggedFile(qscript.dbsnp,"dbsnp,VCF,known=true,training=false,truth=false,prior=2.0")
      this.use_annotation ++= List( "FS", "ReadPosRankSum", "MQRankSum")
      this.maxGaussians = 4
    }
    this.num_threads = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outRecal + ".vqsr"
    this.jobName = queueLogDir + outRecal + ".vqsr"
  }

  case class applyvqsr(inVcf: File, inRecal: File, inTranches: File, outVcf: File, recalmode: org.broadinstitute.sting.gatk.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode) extends ApplyRecalibration with CommandLineGATKArgs {
    this.input = List(inVcf)
    this.recal_file = inRecal
    this.tranches_file = inTranches
    this.out = outVcf
    this.ts_filter_level = 99.9
    this.mode = recalmode
    this.num_threads = 8
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
    this.num_threads = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".combine"
    this.jobName = queueLogDir + outVcf + ".combine"
  }

  case class annotate(inVcf: File, outVcf: File) extends VariantAnnotator with CommandLineGATKArgs {
    this.variant = inVcf
    this.out = outVcf
    this.dbsnp = qscript.dbsnp
    //this parameter seems to have been removed
    //this.requireStrictAlleleMatch = true
    this.group = List("StandardAnnotation")
    //this.useAllAnnotations = true //requires snpEffFile
    this.num_threads = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".annotate"
    this.jobName = queueLogDir + outVcf + ".annotate"
  }

 /****************************************************************************
   * Classes (non-GATK programs)
   ****************************************************************************/

 /*
  case class analyzeCovariates (inRecalFile: File, outPath: File) extends AnalyzeCovariates {
    this.recal_file = inRecalFile
    this.output_dir = outPath.toString
    this.analysisName = queueLogDir + inRecalFile + ".analyze_covariates"
    this.jobName = queueLogDir + inRecalFile + ".analyze_covariates"
  }
*/

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

}

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

import collection.JavaConversions._

import org.apache.commons.io.FilenameUtils

import org.broadinstitute.gatk.engine.arguments.ValidationExclusion
import org.broadinstitute.gatk.engine.filters.ReassignOneMappingQualityFilter
import org.broadinstitute.gatk.engine.phonehome.GATKRunReport
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.function.ListWriterFunction
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.tools.walkers.coverage.DoCOutputType
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypingOutputMode
import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model
import org.broadinstitute.gatk.tools.walkers.genotyper.OutputMode
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.tools.walkers.indels.IndelRealigner.ConsensusDeterminationModel
import org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection //.Mode
import org.broadinstitute.gatk.utils.baq.BAQ.CalculationMode
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType

import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.SAMFileReader
import htsjdk.variant.variantcontext.VariantContext

//http://gatkforums.broadinstitute.org/discussion/1315/frequently-asked-questions-about-scala
//http://www.broadinstitute.org/gatk/guide/tagged?tag=intellij

class Bams2Vcf extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the Bams2Vcf.
  // 'qscript' is now the same as 'Bams2Vcf.this'
  qscript =>

  /****************************************************************************
   * Required Parameters.  All initialized to empty values.
   ****************************************************************************/

  //Required input is one or more of input (-I), extra_recalibrated_bams (-i), and/or existing_gvcfs
  @Input(fullName="input", shortName="I", doc="input aligned, unrecalibrated BAM file(s) with readgroups merged by library, duplicates marked, and de-duped library bam files merged by sample", required=false)
  var input: List[File] = Nil

  @Input(fullName="extra_recalibrated_bams", shortName="i", doc="input realigned recalibrated sample bam file(s) to use as additional (or only) input to Unified Genotyper and subsequent analyses", required=false)
  var extra_recalibrated_bams: List[File] = Nil

  @Input(fullName="rnaseq_unsplit_bams", shortName="r", doc="input unsplit rnaseq sample bam file(s) with readgroups previously added, which need SplitNCigarReads before downstream processing and analyses (ignored if rnaseq is false)", required=false)
  var rnaseq_unsplit_bams: List[File] = Nil

  @Input(fullName="existing_gvcfs", shortName="g", doc="input existing gVCF file(s) to use for joint genotyping with new gVCFs produced by HaplotypeCaller", required=false)
  var existing_gvcfs: List[File] = Nil

  @Argument(fullName="vcfPrefix", shortName="V", doc="prefix (path/name without extensions) of VCF files output by HaplotypeCaller/UnifiedGenotyper and subsequent analyses", required=true)
  var vcf_prefix: String = ""

  /****************************************************************************
   * Optional Parameters. Some initialized to defaults
   ****************************************************************************/

  @Input(fullName="realigned_recalibrated_sample_rename_mapping_file", shortName="sample_rename_mapping_file", doc="sample_rename_mapping_file", required=false)
  var realigned_recalibrated_sample_rename_mapping_file: List[File] = Nil

  @Input(fullName="reference", shortName="R", doc="reference fasta file (default: /share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/human_g1k_v37_decoy.fasta)", required=false)
  var reference: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/human_g1k_v37_decoy.fasta"

  @Input(fullName="dbsnp", shortName="D", doc="dbsnp VCF file (default: /share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/dbsnp_138.b37.vcf)", required=false)
  var dbsnp: File = "/share/apps/myourshaw/resources/gatk_resource_bundle/current/b37/dbsnp_138.b37.vcf"

  //@Input(fullName="bam_interval_string", shortName="bam_interval_string", doc="the -L interval string to be used by TableRecalibration - output recalibrated bams at interval only", required=false)
  //var bam_interval_string: String = ""

  //@Input(fullName="bam_interval_file", shortName="bam_interval_file", doc="an intervals file to be used by quality score recalibration - output variants at intervals only", required=false)
  //var bam_interval_file: File = _

  @Input(fullName="doc_interval_file", shortName="doc_interval_file", doc="an intervals file to be used by DepthOfCoverage - count reads at intervals only (default: /share/apps/myourshaw/resources/intervals/RefGeneIntervals/RefGene.cds.interval_list)", required=false)
  var doc_interval_file: File = "/share/apps/myourshaw/resources/intervals/RefGeneIntervals/RefGene.cds.interval_list"

  @Input(fullName="doc_gene_list", shortName="doc_gene_list", doc="calculate the depth of coverage statistics over this list of genes (default: /share/apps/myourshaw/resources/intervals/RefGeneIntervals/RefGene.gatk.txt)", required=false)
  var doc_gene_list: File = "/share/apps/myourshaw/resources/intervals/RefGeneIntervals/RefGene.gatk.txt"

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

  @Argument(fullName="rnaseq", shortName="rnaseq", doc="input bam files are from mapped RNASeq reads (default: false)", required=false)
  var rnaseq: Boolean = false

  @Argument(fullName="available_nodes", shortName="nodes", doc="number of cluster nodes for scatter/gather (default: 80)", required=false)
  var available_nodes: Int = 80

  @Argument(fullName="scatter_count_max", shortName="scatter_count_max", doc="target maximum number of number of parallel jobs for scatter/gather (default: 0)", required=false)
  //if scatter_count_max argument is < 0, use -1 (no scatter)
  //if scatter_count_max argument is = 0 use max(nContigs, available_nodes)
  //if scatter_count_max argument is > 0 use scatter_count_max
  var scatter_count_max: Int = 0

  /****************************************************************************
   * Hidden Parameters
   ****************************************************************************/
  @Hidden
  @Argument(fullName="scatter_gather", shortName="sg", doc="how many ways to scatter/gather", required=false)
  var nContigs: Int = -1


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

  assert(!(input.isEmpty && extra_recalibrated_bams.isEmpty && rnaseq_unsplit_bams.isEmpty()), "No bam files in input")

   //count number of samples
   val allBamFiles = input ::: extra_recalibrated_bams ::: rnaseq_unsplit_bams
   sampleCount = mapSampleFiles(allBamFiles).keys.size
   //allBamFiles = input ::: extra_recalibrated_bams
   //// map each sample ID to bam file(s)
   //val mapAllBamFiles: Map[String, List[File]] = mapSampleFiles(allBamFiles)
   //sampleCount = mapAllBamFiles.keys.size

   // the number of contigs in the first bam file in the set (may be used to set scatterCount)
   if (nContigs < 0) {
     nContigs = QScriptUtils.getNumberOfContigs(allBamFiles.head)
   }

   //if scatter_count_max argument is < 0, use -1 (no scatter)
   //if scatter_count_max argument is = 0 use max(nContigs, available_nodes)
   //if scatter_count_max argument is > 0 use scatter_count_max
   scatter_count_max = if (scatter_count_max > 0) scatter_count_max else if (scatter_count_max == 0 && nContigs > 1 && available_nodes > 1) nContigs.max(available_nodes) else -1

  if(rnaseq) {
    haplotype_caller = true
    unified_genotyper = false
    if (!rnaseq_unsplit_bams.isEmpty()){
      //SplitNCigarReads
      // uniquify list of rnaseq unsplit bam files
      val unsplitBams = rnaseq_unsplit_bams.toSet

      //try not to create excessive scatter-gathering
      val unsplitBams_scatter_count = scatter_count_max / unsplitBams.size

      //file names for output derived from input bam files
      for (unsplitbam <- unsplitBams) {

        val splitBam = swapExt(unsplitbam, ".bam", ".split.bam")

        add(splitNCigarReads(List(unsplitbam), splitBam))

        input :+ splitBam
      }
    }
  }

  // uniqueify input list of bam files that are already recalibrated
  var recalibratedBams = extra_recalibrated_bams.toSet

  //uniqueify input list of single-sample or combined gVCF files already produced by HaplotypeCsller
  var gvcf_hcs = existing_gvcfs.toSet

   // realign and recalibrate raw bam files, if any
  if (!input.isEmpty) {

    // uniqueify input list of bam files
    val unrecalibratedBams = input.toSet

    //try not to create excessive scatter-gathering
    val unrecalibratedBams_scatter_count = scatter_count_max / unrecalibratedBams.size

    // set the model for the Indel Realigner
    cleanModelEnum = getIndelCleaningModel()

    // if this is a 'knowns only' indel realignment run, do it only once for all samples.
    val globalIntervals = new File(outputDir + qscript.projectName + ".intervals")
    if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY)
      add(realignerTargetCreator(null, globalIntervals, unrecalibratedBams_scatter_count))

    // map each sample ID to bam file(s)
    val sampleBAMFiles: Map[String, List[File]] = mapSampleFiles(unrecalibratedBams.toList)

    //file names for output derived from input bam files
    for ((sample, bamList) <- sampleBAMFiles) {

      // BAM files generated by the pipeline
      //new file name if a sample is in multiple bam files
      val bam = if (bamList.length == 1) bamList(0) else new File(qscript.projectName + "." + sample + ".bam")
      val realignedBam = swapExt(bam, ".bam", ".realigned.bam")
      val dedupedBam = swapExt(realignedBam, ".bam", ".dedup.bam")
      val recalibratedBam   = if (dedup_cleaned_bams) swapExt(bam, ".bam", ".realigned.dedup.recalibrated.bam") else swapExt(bam, ".bam", ".realigned.recalibrated.bam")

      // Accessory files
      val targetIntervals = if (cleaningModel == ConsensusDeterminationModel.KNOWNS_ONLY) {globalIntervals} else {swapExt(bam, ".bam", ".intervals")}
      val dedupMetricsFile = swapExt(bam, ".bam", ".dedup_metrics")
      val recal_dataFile = swapExt(bam, ".bam", ".recal_data.first_pass.table")
      val recal_dataFile2 = swapExt(bam, ".bam", ".recal_data.second_pass.table")
      val recal_plots = swapExt(bam, ".bam", ".recal.plots")
      val preValidateLog = swapExt(bam, ".bam", ".pre.validation")
      val postValidateLog = swapExt(bam, ".bam", ".post.validation")

      // Validation is an optional step for the BAM file generated after
      // alignment and the final bam file of the pipeline.
      /*
        rnaseq files will WARN all records with
        NM tag (nucleotide differences) is missing
        and some records will ERROR
        No real operator (M|I|D|N) in CIGAR
      */
      if (!noValidation) {
        if (rnaseq) {
          for (sampleFile <- bamList)
            add(validate(sampleFile, preValidateLog, List("MISSING_TAG_NM", "INVALID_CIGAR")),
              validate(recalibratedBam, postValidateLog, List("MISSING_TAG_NM", "INVALID_CIGAR")))
        } else {
          for (sampleFile <- bamList)
            add(validate(sampleFile, preValidateLog),
              validate(recalibratedBam, postValidateLog))
        }
      }

      //IndelRealigner
      if (cleaningModel != ConsensusDeterminationModel.KNOWNS_ONLY)
        add(realignerTargetCreator(bamList, targetIntervals, unrecalibratedBams_scatter_count))
      add(indelRealigner(bamList, targetIntervals, realignedBam, unrecalibratedBams_scatter_count))

      //BaseRecalibrator optionally dedup and then generate a first pass recalibration table file
      if (dedup_cleaned_bams){
        add(dedup(realignedBam, dedupedBam, dedupMetricsFile))
        add(baseRecalibrator(List(dedupedBam), recal_dataFile, unrecalibratedBams_scatter_count))
      }
      else{
        add(baseRecalibrator(List(realignedBam), recal_dataFile, unrecalibratedBams_scatter_count))
      }
      //AnalyzeCovariates by generating a second pass recalibration table file
      add(baseRecalibrator(List(realignedBam), recal_dataFile2, unrecalibratedBams_scatter_count, Some(recal_dataFile)))
      add(analyzeCovariates(recal_dataFile, recal_dataFile2, recal_plots))

      // create a recalibrated bam file
      add(createRecalibrated(List(realignedBam), recal_dataFile, recalibratedBam))
      //add(cov(recalibratedBam, postRecalFile),
      //analyzeCovariates(preRecalFile, preOutPath),
      //analyzeCovariates(postRecalFile, postOutPath))

      recalibratedBams += recalibratedBam
    }
  }

  // union of new and extra recalibrated bams
  var recalibratedSampleBams: List[File] = List()


  //ReduceReads obsolete as of GATK 3.0
  /*
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
  */

  recalibratedSampleBams = recalibratedBams.toList.sortWith(_.toString < _.toString)

  val recalibratedSampleBams_scatter_count = scatter_count_max / recalibratedSampleBams.size

  // variant genotyping and annotation on recalibrated input bams + previously recalibrated extra bams
  if (!recalibratedSampleBams.isEmpty){

    //depth of coverage
    val doc_output_prefix =  qscript.vcf_prefix + ".coverage"

    //file names for HaplotypeCaller
    //gVCFs (must end with .vcf)
    var new_existing_gvcfs: List[File] = List()
    val raw_gvcf_hc = qscript.vcf_prefix + ".hc.raw.gvcf.vcf"
    //raw
    val raw_vcf_hc = qscript.vcf_prefix + ".hc.raw.vcf"
    //eval raw
    val raw_vcf_eval_hc = swapExt(raw_vcf_hc, ".vcf",".VariantEval")
    //val hc_metrics = swapExt(raw_vcf_hc, ".vcf", ".metrics")
    //recalibrate snvs
    val raw_snvs_vcf_hc = vcf_prefix + ".hc.raw.snvs.vcf"
    val vqsr_snvs_recal_hc = raw_snvs_vcf_hc + ".recal"
    val vqsr_snvs_tranches_hc = raw_snvs_vcf_hc + ".tranches"
    val vqsr_snvs_plots_hc = raw_snvs_vcf_hc + "plots.R"
    val recalibrated_snvs_vcf_hc = qscript.vcf_prefix + ".hc.recalibrated.snvs.vcf"
    val filtered_snvs_vcf_hc = qscript.vcf_prefix + ".hc.filtered.snvs.vcf" //RNASeq
    //recalibrate indels
    val raw_indels_vcf_hc = qscript.vcf_prefix  + ".hc.raw.indels.vcf"
    val vqsr_indels_recal_hc = raw_indels_vcf_hc + ".recal"
    val vqsr_indels_tranches_hc = raw_indels_vcf_hc + ".tranches"
    val vqsr_indels_plots_hc = raw_indels_vcf_hc + "plots.R"
    val recalibrated_indels_vcf_hc = qscript.vcf_prefix + ".hc.recalibrated.indels.vcf"
    val filtered_indels_vcf_hc = qscript.vcf_prefix + ".hc.filtered.indels.vcf" //RNASeq
    //combine snvs & indels
    val analysis_ready_vcf_hc = qscript.vcf_prefix + ".hc.analysis_ready.vcf"
    //eval
    val analysis_ready_vcf_eval_hc = swapExt(analysis_ready_vcf_hc, ".vcf",".VariantEval")
    //annotate
    val annotated_vcf_hc = swapExt(analysis_ready_vcf_hc, ".vcf",".annotated.vcf")

    //file names for UnifiedGenotyper
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

    //val indel_filter_expression = if (sampleCount >= 10) List("QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0") else List("QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0")
    val rnaseq_filter_names = List("RNASeq_QDFS")
    val rnaseq_filter_expressions = List("QD < 2.0 || FS > 30.0")

    //HaplotypeCaller joint discovery
    //joint discovery workflow is deprecated; use incremental joint variant discovery instead
    /*
    if(haplotype_caller){
      add(haplotypeCaller(recalibratedSampleBams, raw_vcf_hc, recalibratedSampleBams_scatter_count_max))
      //split snvs and indels into two files
      add(selectVariants(raw_vcf_hc, raw_snvs_vcf_hc, List(VariantContext.Type.SNP))
        ,selectVariants(raw_vcf_hc, raw_indels_vcf_hc, List(VariantContext.Type.INDEL)))
      //VQSR
      add(vqsr(raw_snvs_vcf_hc, vqsr_snvs_recal_hc, vqsr_snvs_tranches_hc, vqsr_snvs_plots_hc, VariantRecalibratorArgumentCollection.Mode.SNP, false)
        ,applyvqsr(raw_snvs_vcf_hc, vqsr_snvs_recal_hc, vqsr_snvs_tranches_hc, recalibrated_snvs_vcf_hc, VariantRecalibratorArgumentCollection.Mode.SNP))
      add(vqsr(raw_indels_vcf_hc, vqsr_indels_recal_hc, vqsr_indels_tranches_hc, vqsr_indels_plots_hc, VariantRecalibratorArgumentCollection.Mode.INDEL, false)
        ,applyvqsr(raw_indels_vcf_hc, vqsr_indels_recal_hc, vqsr_indels_tranches_hc, recalibrated_indels_vcf_hc, VariantRecalibratorArgumentCollection.Mode.INDEL))
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
    */
    //HaplotypeCaller incremental joint variant discovery
    if(haplotype_caller){

      //create single-sample gVCFs
      for (recalibratedSampleBam <- recalibratedSampleBams){
        val gvcf_hc = recalibratedSampleBam + ".hc.gvcf.vcf"
        gvcf_hcs += gvcf_hc
        add(haplotypeCaller(List(recalibratedSampleBam), gvcf_hc, recalibratedSampleBams_scatter_count))
      }

      //hierarchically merge new and existing gVCFs into a single gVCF for input to GenotypeGVCFs
      new_existing_gvcfs = gvcf_hcs.toList.sortWith(_.toString < _.toString)
      add(combineGVCFs(new_existing_gvcfs, raw_gvcf_hc, scatter_count_max))

      //joint genotyping to produce raw SNP and indel VCFs for input to VQSR
      add(genotypeGVCFs(List(raw_gvcf_hc), raw_vcf_hc, scatter_count_max))

        //split snvs and indels into two files
        add(selectVariants(raw_vcf_hc, raw_snvs_vcf_hc, List(VariantContext.Type.SNP))
          ,selectVariants(raw_vcf_hc, raw_indels_vcf_hc, List(VariantContext.Type.INDEL)))

        if (rnaseq) {
          //hard filters
          add(filter(raw_snvs_vcf_hc, filtered_snvs_vcf_hc, rnaseq_filter_names, rnaseq_filter_expressions))
          add(filter(raw_indels_vcf_hc, filtered_indels_vcf_hc, rnaseq_filter_names, rnaseq_filter_expressions))
          add(combine(List(filtered_snvs_vcf_hc, filtered_indels_vcf_hc), analysis_ready_vcf_hc))
        } else {
        //VQSR
        add(vqsr(raw_snvs_vcf_hc, vqsr_snvs_recal_hc, vqsr_snvs_tranches_hc, vqsr_snvs_plots_hc, VariantRecalibratorArgumentCollection.Mode.SNP, false)
          ,applyvqsr(raw_snvs_vcf_hc, vqsr_snvs_recal_hc, vqsr_snvs_tranches_hc, recalibrated_snvs_vcf_hc, VariantRecalibratorArgumentCollection.Mode.SNP))
        add(vqsr(raw_indels_vcf_hc, vqsr_indels_recal_hc, vqsr_indels_tranches_hc, vqsr_indels_plots_hc, VariantRecalibratorArgumentCollection.Mode.INDEL, false)
          ,applyvqsr(raw_indels_vcf_hc, vqsr_indels_recal_hc, vqsr_indels_tranches_hc, recalibrated_indels_vcf_hc, VariantRecalibratorArgumentCollection.Mode.INDEL))

          //combine snvs and indels
          add(combine(List(recalibrated_snvs_vcf_hc, recalibrated_indels_vcf_hc), analysis_ready_vcf_hc))
      }

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
      add(unifiedGenotyper(recalibratedSampleBams, raw_vcf_ug, scatter_count_max))

      //split snvs and indels into two files
      add(selectVariants(raw_vcf_ug, raw_snvs_vcf_ug, List(VariantContext.Type.SNP))
        ,selectVariants(raw_vcf_ug, raw_indels_vcf_ug, List(VariantContext.Type.INDEL)))

      //VQSR
      add(vqsr(raw_snvs_vcf_ug, vqsr_snvs_recal_ug, vqsr_snvs_tranches_ug, vqsr_snvs_plots_ug, VariantRecalibratorArgumentCollection.Mode.SNP, true)
        ,applyvqsr(raw_snvs_vcf_ug, vqsr_snvs_recal_ug, vqsr_snvs_tranches_ug, recalibrated_snvs_vcf_ug, VariantRecalibratorArgumentCollection.Mode.SNP))
      add(vqsr(raw_indels_vcf_ug, vqsr_indels_recal_ug, vqsr_indels_tranches_ug, vqsr_indels_plots_ug, VariantRecalibratorArgumentCollection.Mode.INDEL, true)
        ,applyvqsr(raw_indels_vcf_ug, vqsr_indels_recal_ug, vqsr_indels_tranches_ug, recalibrated_indels_vcf_ug, VariantRecalibratorArgumentCollection.Mode.INDEL))
      /*
      if (!is_whole_genome){
        add(filter(raw_indels_vcf_ug, filtered_indels_vcf_ug, indel_filter_expression, List("GATKStandard")))
        add(combine(List(recalibrated_snvs_vcf_ug, filtered_indels_vcf_ug), analysis_ready_vcf_ug))
      }
      else{
        add(vqsr(raw_indels_vcf_ug, vqsr_indel_recal_ug, indel_tranches_ug, vqsr_indel_plots_ug, VariantRecalibratorArgumentCollection.Mode.INDEL)
          ,applyvqsr(raw_indels_vcf_ug, vqsr_indel_recal_ug, indel_tranches_ug, recalibrated_indel_vcf_ug, VariantRecalibratorArgumentCollection.Mode.INDEL))
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

    if(depth_of_coverage){
      //scatter-gather may only work with depth by locus
      add(depth(recalibratedSampleBams, doc_output_prefix))
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

  case class realignerTargetCreator (inBams: List[File], outIntervals: File, scatter_count: Int) extends RealignerTargetCreator with CommandLineGATKArgs {
    if (cleanModelEnum != ConsensusDeterminationModel.KNOWNS_ONLY)
      this.input_file = inBams
    this.out = outIntervals
    //0.0 is default
    //this.mismatchFraction = 0.0
    this.known =List(qscript.mills, qscript.onekg_indels)
    this.scatterCount = scatter_count
    this.num_threads = 8 //4 maybe 6 or even 8?
    this.isIntermediate = false
    this.analysisName = queueLogDir + outIntervals + ".target"
    this.jobName = queueLogDir + outIntervals + ".target"
  }

  case class indelRealigner (inBams: List[File], tIntervals: File, outBam: File, scatter_count: Int) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file = inBams
    this.targetIntervals = tIntervals
    this.out = outBam
    this.known = List(qscript.mills, qscript.onekg_indels)
    //USE_READS is default
    this.consensusDeterminationModel = ConsensusDeterminationModel.USE_READS
//    this.compress = 0
//    this.noPGTag = qscript.testMode;
    this.scatterCount = scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".clean"
    this.jobName = queueLogDir + outBam + ".clean"
  }

  case class baseRecalibrator(inBams: List[File], recal_data: File, scatter_count: Int, BQSR_recal_data: Option[String] = None) extends BaseRecalibrator with CommandLineGATKArgs {
    this.input_file = inBams
    //None for first pass, first pass recal_data for second pass
    if (! BQSR_recal_data.isEmpty) this.BQSR = BQSR_recal_data.get
    this.knownSites ++= List(qscript.dbsnp, qscript.mills, qscript.onekg_indels)
    this.out = recal_data
    //this parameter seems to have been dropped
    //omit plot_pdf_file if using Queue scatter/gather
    //this.plot_pdf_file = recal_data+".plot.pdf"
    //not sure if this is needed, it may relate to an older version of Novoalign
    this.read_filter = List("BadCigar")
    //this should work but has failed in gather combine tables
    //BQSRGatherer: ... java.lang.IllegalArgumentException: Table1 2,3 not equal to 1,3
    //this.scatterCount = scatter_count
    //http://gatkforums.broadinstitute.org/discussion/3265/bqsrgatherer-exception
    //try with 20
    this.scatterCount = scatter_count.max(20)
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
    this.num_cpu_threads_per_data_thread = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".createRecalibrated"
    this.jobName = queueLogDir + outBam + ".createRecalibrated"
  }

  case class analyzeCovariates(beforeTable: File, afterTable: File, plotsPdf: File) extends AnalyzeCovariates with CommandLineGATKArgs {
    this.before = beforeTable
    this.after = afterTable
    this.plots = plotsPdf
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
    this.isIntermediate = false
    this.analysisName = queueLogDir + outPrefix + ".depth"
    this.jobName = queueLogDir + outPrefix + ".depth"
  }

  case class haplotypeCaller (inBams: List[File], outVcf: File, scatter_count: Int) extends HaplotypeCaller with CommandLineGATKArgsHiMem {
    //the default memoryLimit of 6 may be OK for ReduceReads BAMs more may be required for uncompressed BAMs
    //try default memoryLimit with single-sample workflow
    this.input_file = inBams
    this.out = outVcf
    this.intervals = List(qscript.ug_interval_file)
    this.dbsnp = qscript.dbsnp
    this.read_filter = List("BadCigar")
    this.genotyping_mode = GenotypingOutputMode.DISCOVERY
    if (rnaseq) {
      this.standard_min_confidence_threshold_for_calling = 20.0
      this.standard_min_confidence_threshold_for_emitting = 20.0
    } else {
      this.standard_min_confidence_threshold_for_calling = 30.0
      this.standard_min_confidence_threshold_for_emitting = 10.0
    }
    this.downsample_to_coverage = 250
    this.group = List("Standard")
    //next three parameters are new to GATK 3.0
    this.emitRefConfidence = ReferenceConfidenceMode.GVCF
    this.variant_index_type = GATKVCFIndexType.LINEAR
    this.variant_index_parameter = 128000
    if (rnaseq) {
      this.recoverDanglingHeads = true
      this.dontUseSoftClippedBases = true
    }
    //##### ERROR MESSAGE: Invalid command line: Argument baq has a bad value: Walker cannot accept BAQ'd base qualities, and yet BAQ mode CALCULATE_AS_NECESSARY was requested.
    //this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    //##### ERROR MESSAGE: Invalid command line: Argument nt has a bad value: The analysis HaplotypeCaller currently does not support parallel execution with nt.  Please run your analysis without the nt option.
    //this.num_threads = 8
    //this.jobNativeArgs  = Seq("-w n -l excl=true -l mem_free=20G")
    //num_cpu_threads_per_data_thread should be <= number of cores
    this.num_cpu_threads_per_data_thread = 8
    this.scatterCount = scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".hc"
    this.jobName = queueLogDir + outVcf + ".hc"
  }

  case class combineGVCFs(inGvcfs: List[File], outGvcf: File, scatter_count: Int) extends CombineGVCFs with CommandLineGATKArgs {
    this.variant = inGvcfs
    this.out = outGvcf
    //best parallelism?
    //nt not supported
    //this.num_threads = 8
    //nct not supported
    //this.num_cpu_threads_per_data_thread = 8
    this.scatterCount = scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outGvcf + ".combine"
    this.jobName = queueLogDir + outGvcf + ".combine"
  }

  case class genotypeGVCFs (inGvcfs: List[File], outVcf: File, scatter_count: Int) extends GenotypeGVCFs with CommandLineGATKArgsHiMem {
    this.variant = inGvcfs
    this.out = outVcf
    this.dbsnp = qscript.dbsnp
    //defaults:
    this.annotation = List("InbreedingCoeff", "FisherStrand", "QualByDepth", "ChromosomeCounts")
    this.num_threads = 8
    this.jobNativeArgs  = Seq("-w n -l excl=true -l mem_free=20G")
    //nct not supported
    //this.num_cpu_threads_per_data_thread = 8
    this.scatterCount = scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".combine"
    this.jobName = queueLogDir + outVcf + ".combine"
  }

  //case class unifiedGenotyper (inBams: List[File], outVcf: File, metricsFile: File) extends UnifiedGenotyper with CommandLineGATKArgs {
  case class unifiedGenotyper (inBams: List[File], outVcf: File, scatter_count: Int) extends UnifiedGenotyper with CommandLineGATKArgs {
    // If you are running this on a compute farm, make sure that the Sting/shell
    // folder is in your path to use mergeText.sh and splitIntervals.sh.
    this.input_file = inBams
    this.out = outVcf
    //this.metrics_file = metricsFile
    this.intervals = List(qscript.ug_interval_file)
    this.dbsnp = qscript.dbsnp
    this.read_filter = List("BadCigar")
    this.output_mode = OutputMode.EMIT_VARIANTS_ONLY
    this.genotype_likelihoods_model = Model.BOTH
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
    this.scatterCount = scatter_count
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".ug"
    this.jobName = queueLogDir + outVcf + ".ug"
  }

  case class variantEval (inVcf: File, outEval: File) extends VariantEval with CommandLineGATKArgs {
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

  case class selectVariants(inVcf: File, outVcf: File, seltype: List[VariantContext.Type]) extends SelectVariants with CommandLineGATKArgs {
    this.variant = inVcf
    this.out = outVcf
    this.selectTypeToInclude = seltype
    this.num_threads = 8
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".select"
    this.jobName = queueLogDir + outVcf + ".select"
  }

  case class vqsr(inVcf: File, outRecal: File, outTranches: File, outPlots: File, recalmode: VariantRecalibratorArgumentCollection.Mode, input_from_ug: Boolean) extends VariantRecalibrator with CommandLineGATKArgs {
    this.input :+= inVcf
    this.recal_file = outRecal
    this.tranches_file = outTranches
    this.rscript_file = outPlots
    this.mode = recalmode
    //percentBad obsolete as of GATK 2.7
    //this.percentBad = if(sampleCount >= 30) 0.01 else 0.05
    this.minNumBadVariants = 1000
    if (is_whole_genome) {
      this.use_annotation ++= List("DP")
    }
    if (sampleCount >= 10){
      this.use_annotation ++= List("InbreedingCoeff")
    }
    if (mode == VariantRecalibratorArgumentCollection.Mode.SNP) {
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
    else if (mode == VariantRecalibratorArgumentCollection.Mode.INDEL) {
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

  case class applyvqsr(inVcf: File, inRecal: File, inTranches: File, outVcf: File, recalmode: VariantRecalibratorArgumentCollection.Mode) extends ApplyRecalibration with CommandLineGATKArgs {
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
    if (rnaseq) {
      this.window = 35
      this.cluster = 3
    }
    this.isIntermediate = false
    this.analysisName = queueLogDir + outVcf + ".filter"
    this.jobName = queueLogDir + outVcf + ".filter"
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
    * Classes RNASeq
    ****************************************************************************/
  // http://www.broadinstitute.org/gatk/guide/article?id=3891

  //may require
  //case class SplitNCigarReads (inBam: File, outBam: File) extends SplitNCigarReads with CommandLineGATKArgsHiMem {
  case class splitNCigarReads (inBam: List[File], outBam: File) extends SplitNCigarReads with CommandLineGATKArgs {
    this.I = inBam
    this.out = outBam
    //the default DMQ is 60
    this.read_filter = List("ReassignMappingQuality -DMQ 60")
    this.unsafe = ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS
    this.isIntermediate = false
    this.analysisName = queueLogDir + outBam + ".splitNCigar"
    this.jobName = queueLogDir + outBam + ".splitNCigar"
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

  case class validate (inBam: File, outLog: File, ignore: List[String] = Nil) extends ValidateSamFile with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outLog
    this.REFERENCE_SEQUENCE = qscript.reference
    if (!ignore.isEmpty()){
      this.IGNORE = ignore
    }
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

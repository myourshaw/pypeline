package org.broadinstitute.sting.queue.qscripts.NL

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.gatk.phonehome.GATKRunReport

import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.commandline.Hidden

import collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import org.apache.commons.io.FilenameUtils

/**
 * An example building on the intro ExampleCountReads.scala.
 * Runs an INCOMPLETE version of the UnifiedGenotyper with VariantEval and optional VariantFiltration.
 */
class ExomeUnifiedGenotyper extends QScript {
  // Create an alias 'qscript' to be able to access variables
  // in the ExampleUnifiedGenotyper.
  // 'qscript' is now the same as 'ExampleUnifiedGenotyper.this'
  qscript =>


  // Required arguments.  All initialized to empty values.

  @Input(doc="input BAM file - or list of BAM files (can be specified mulitple times)", fullName="input", shortName="i", required=true)
  var input: List[File] = _

  @Input(doc="The reference file for the bam files.", fullName="reference", shortName="R", required=true)
  var reference: File = _ // _ is scala shorthand for null

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: File = _

  @Input(doc="hapmap ROD to use (must be in VCF format)", fullName="hapmap", shortName="H", required=true)
  var hapmap: File = _

  @Input(doc="omni ROD to use (must be in VCF format)", fullName="omni", shortName="omni", required=true)
  var omni: File = _

  @Input(doc="1000 Genomes project consensus ROD to use (must be in VCF format)", fullName="kg", shortName="kg", required=true)
  var kg: File = _
  // The following arguments are all optional.

  @Input(doc="the project name determines the final output (vcf file) base name. Example NA12878 yields NA12878.vcf", fullName="project", shortName="p", required=false)
  var projectName: String = "project"

  @Argument(shortName="outputDir", doc="output directory", required=false)
  var outputDir: String = "./"

  @Argument(shortName="noBAQ", doc="turns off BAQ calculation", required=false)
  var noBAQ: Boolean = false

  @Argument(shortName="noIndels", doc="do not call indels with the Unified Genotyper", required=false)
  var noIndels: Boolean = false

  @Input(doc="An optional file with a list of intervals to proccess.", shortName="L", required=false)
  var intervals: File = _

  @Input(doc="An optional file with a list of intervals to evaluate.", shortName="eval-intervals", required=false)
  var eval_intervals: File = _

  @Argument(shortName="ST", doc="list of stratification modules", required=false)
  var ST: List[String] = Nil

  @Argument(doc="A optional list of filter names.", shortName="filter", required=false)
  var filterNames: List[String] = Nil // Nil is an empty List, versus null which means a non-existent List.

  @Argument(doc="An optional list of filter expressions.", shortName="filterExpression", required=false)
  var filterExpressions: List[String] = Nil

  @Argument(doc="Ti/tv target (used only for plotting).", shortName="titv", required=false)
  var titv: Double = 3.0

  @Argument(shortName="mbq", doc="The minimum Phred-Scaled quality score threshold to be considered a good base.", required=false)
  var minimumBaseQuality: Int = -1

  @Argument(shortName="deletions", doc="Maximum deletion fraction allowed at a site to call a genotype.", required=false)
  var deletions: Double = -1

  @Argument(shortName="alleles", doc="The set of alleles at which to genotype (optional).", required=false)
  var alleles: File = _


  /****************************************************************************
  * Hidden Parameters
  *************************************************************************** */
  @Hidden
  @Input(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = -1

  //val training_1000G = "/home/kmsquire/Data/resources/ALL.wgs.phase1.projectConsensus.snps.sites.ALL_PASS.vcf.gz"
  //val badSites_1000G = "/humgen/1kg/processing/official_release/phase1/projectConsensus/phase1.wgs.projectConsensus.v2b.recal.terrible.vcf"
  //val projectConsensus_1000G = "/home/kmsquire/Data/resources/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.gz"

  val queueLogDir = ".qlog/"

  class Target(
          val baseName: String,
          val reference: File,
          val dbsnpFile: File,
          val hapmapFile: File,
          val omniFile: File,
          val kgFile: File,
          val bamList: List[File],
          val sampleList: List[String],
          val intervals: File,
          val trancheTarget: Double,
          val nContigs: Int,
          val nSamples: Int) {
    val isExome = false
    val name = FilenameUtils.concat(qscript.outputDir, baseName)
    val clusterFile = new File(name + ".clusters")
    val rawVCF = new File(name + ".raw.vcf")
    val rawIndelVCF = new File(name + ".raw.indel.vcf")
    val filteredIndelVCF = new File(name + ".filtered.indel.vcf")
    val recalibratedVCF = new File(name + ".recalibrated.vcf")
    val tranchesFile = new File(name + ".tranches")
    val vqsrRscript = name + ".vqsr.r"
    val recalFile = new File(name + ".tranches.recal")
    val evalFile = new File(name + ".snp.eval")
    val evalIndelFile = new File(name + ".indel.eval")
  }

  def getSampleNames(bamFiles: List[File]): List[String] = {
    var sampleList = List.empty[String]
    for (bam <- bamFiles) {
      val samReader = new SAMFileReader(bam)
      val header = samReader.getFileHeader
      val readGroups = header.getReadGroups

      // Fill out the sample list with the samples in this bam
      for (rg <- readGroups) {
        val sample = rg.getSample
        if(!sampleList.contains(sample))
          sampleList +: sample
      }
    }
    return sampleList
  }



  def script = {

    // keep a record of the number of contigs in the first bam file in the list
    val bams = 
      if (input.length == 1)
        QScriptUtils.createListFromFile(input(0))
      else
        input
    if (nContigs < 0)
     nContigs = QScriptUtils.getNumberOfContigs(bams(0))

    val samples = getSampleNames(bams)

    val target = new Target(
      qscript.projectName,
      qscript.reference,
      qscript.dbSNP,
      qscript.hapmap,
      qscript.omni,
      qscript.kg,
      bams,
      samples,
      qscript.intervals,
      qscript.titv,
      nContigs,
      bams.size
    )

    if (!noIndels) add(indelCall(target), indelFilter(target), indelEvaluation(target))
    add(snpCall(target))
    add(VQSR(target))
    add(applyVQSR(target))
    add(snpEvaluation(target))

  }

  // General arguments to GATK walkers
  trait CommandLineGATKArgs extends CommandLineGATK {
    this.reference_sequence = qscript.reference
    this.intervals = if (qscript.intervals == null) Nil else List(qscript.intervals)
    this.logging_level = "INFO"
    this.memoryLimit = 6
    this.phone_home = GATKRunReport.PhoneHomeOption.STANDARD
  }

  def bai(bam: File) = new File(bam + ".bai")

  // 1.) Unified Genotyper Base
  class GenotyperBase (t: Target) extends UnifiedGenotyper with CommandLineGATKArgs {
    this.scatterCount = t.nContigs
    this.nt = 8
    this.dcov = 1000
    this.stand_call_conf = 30.0
    this.stand_emit_conf = 10.0
    this.input_file ++= t.bamList
    this.D = new File(t.dbsnpFile)
  }


  // 1a.) Call SNPs with UG
  case class snpCall (t: Target) extends GenotyperBase(t) {
    if (minimumBaseQuality >= 0)
      this.min_base_quality_score = minimumBaseQuality
    if (qscript.deletions >= 0)
      this.max_deletion_fraction = qscript.deletions
    this.out = t.rawVCF
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.SNP
    this.baq = if (noBAQ) {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF} else {org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.CALCULATE_AS_NECESSARY}
    if (qscript.alleles != null) {
      this.genotyping_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
      this.alleles = qscript.alleles
    }
    this.analysisName = t.name + "_UGs"
    this.jobName =  queueLogDir + t.name + ".snpcall"
  }

  // 1b.) Call Indels with UG
  case class indelCall (t: Target) extends GenotyperBase(t) {
    this.memoryLimit = 6
    this.out = t.rawIndelVCF
    this.glm = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model.INDEL
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.OFF
    if (qscript.alleles != null) {
      this.genotyping_mode = org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES
      this.alleles = qscript.alleles
    }
    this.analysisName = t.name + "_UGi"
    this.jobName =  queueLogDir + t.name + ".indelcall"
  }

  // 2.) Hard Filtering for indels
  case class indelFilter (t: Target) extends VariantFiltration with CommandLineGATKArgs {
    this.memoryLimit = 2
    this.reference_sequence = t.reference
    this.scatterCount = 10
    this.V = t.rawIndelVCF
    this.out = t.filteredIndelVCF
    this.filterName ++= List("QDFilter", "ReadPosFilter", "FSFilter")
    this.filterExpression ++= List("\"QD < 2.0\"", "\"ReadPosRankSum < -20.0\"", "\"FS > 200.0\"")
    //if (t.nSamples >= 10) {
    //    this.filterName ++= List("IndelInbreedingCoeff")
    //    this.filterExpression ++= List("\"InbreedingCoeff < -0.8\"")
    //}
    this.analysisName = t.name + "_VF"
    this.jobName =  queueLogDir + t.name + ".indelfilter"
  }

  // 3.) Variant Quality Score Recalibration - Generate Recalibration table
  case class VQSR(t: Target) extends VariantRecalibrator with CommandLineGATKArgs {
    this.nt = 8
    this.reference_sequence = t.reference
    this.input :+= t.rawVCF
    this.resource :+= new TaggedFile( t.hapmapFile, "training=true,truth=true,prior=15.0" )
    this.resource :+= new TaggedFile( t.omniFile, "training=true,truth=true,prior=12.0" )
    //this.resource :+= new TaggedFile( training_1000G, "training=true,prior=10.0" )
    this.resource :+= new TaggedFile( t.dbsnpFile, "known=true,prior=2.0" )
    //this.resource :+= new TaggedFile( projectConsensus_1000G, "prior=8.0" )
    this.resource :+= new TaggedFile( t.kgFile, "prior=8.0" )
    this.use_annotation ++= List("QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ", "FS")
    if(t.nSamples >= 10) { // InbreedingCoeff is a population-wide statistic that requires at least 10 samples to calculate
        this.use_annotation ++= List("InbreedingCoeff")
    }
    if(!t.isExome) {
        this.use_annotation ++= List("DP")
    } else { // exome specific parameters
        // TODO: make this
        //this.resource :+= new TaggedFile( badSites_1000G, "bad=true,prior=2.0" )
        this.mG = 10
        if(t.nSamples <= 6) { // very few exome samples means very few variants
            this.mG = 4
            this.percentBad = 0.04
        }
    }
    this.tranches_file = t.tranchesFile
    this.recal_file = t.recalFile
    this.allPoly = true
    // Using defaults...?
    this.tranche ++= List("100.0", "99.9", "99.0", "90.0")
    this.rscript_file = t.vqsrRscript
    this.analysisName = t.name + "_VQSR"
    this.jobName = queueLogDir + t.name + ".VQSR"
  }

  // 4.) Apply the recalibration table to the appropriate tranches
  case class applyVQSR (t: Target) extends ApplyRecalibration with CommandLineGATKArgs {
    this.memoryLimit = 6
    this.reference_sequence = t.reference
    this.input :+= t.rawVCF
    this.tranches_file = t.tranchesFile
    this.recal_file = t.recalFile
    this.ts_filter_level = t.trancheTarget
    this.out = t.recalibratedVCF
    this.analysisName = t.name + "_AVQSR"
    this.jobName = queueLogDir + t.name + ".applyVQSR"
  }

  // 5.) Variant Evaluation Base(OPTIONAL)
  class EvalBase(t: Target) extends VariantEval with CommandLineGATKArgs {
    this.memoryLimit = 3
    this.reference_sequence = t.reference
    this.comp :+= new TaggedFile(t.hapmapFile, "hapmap" )
    this.D = new File(t.dbsnpFile)
    this.sample = t.sampleList
  }

  // 5a.) SNP Evaluation (OPTIONAL) based on the cut vcf
  case class snpEvaluation(t: Target) extends EvalBase(t) {
    this.comp :+= new TaggedFile( t.omniFile, "omni" )
    this.eval :+= t.recalibratedVCF
    this.out =  t.evalFile
    this.analysisName = t.name + "_VEs"
    this.jobName = queueLogDir + t.name + ".snp.eval"
    if (qscript.eval_intervals != null)
      this.intervals = List(qscript.eval_intervals)
    if (qscript.ST != Nil)
      this.ST = ST
  }

  // 5b.) Indel Evaluation (OPTIONAL)
  case class indelEvaluation(t: Target) extends EvalBase(t) {
    this.eval :+= t.filteredIndelVCF
    this.evalModule :+= "IndelStatistics"
    this.out =  t.evalIndelFile
    this.analysisName = t.name + "_VEi"
    this.jobName = queueLogDir + queueLogDir + t.name + ".indel.eval"
    if (qscript.eval_intervals != null)
      this.intervals = List(qscript.eval_intervals)
  }


}

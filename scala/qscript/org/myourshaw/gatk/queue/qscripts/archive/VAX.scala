/**
 * Created by IntelliJ IDEA.
 * User: myourshaw
 * Date: 2012-02-18
 * Time: 15:24
 * To change this template use File | Settings | File Templates.
 */

package org.myourshaw.sting.queue.qscripts

import org.broadinstitute.sting.queue.function.CommandLineFunction
import java.io.File
import org.broadinstitute.sting.commandline.Argument
import org.broadinstitute.sting.commandline.Gather
import org.broadinstitute.sting.commandline.Input
import org.broadinstitute.sting.commandline.Output
import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction
import org.broadinstitute.sting.queue.util.ShellUtils
import org.broadinstitute.sting.queue.extensions.gatk.LocusScatterFunction

class VAX extends org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK with ScatterGatherableFunction {
  analysisName = "VAX"
  analysis_type = "VAX"
  scatterClass = classOf[LocusScatterFunction]

  @Input(fullName="input", shortName="i", doc="input VCF file", required=true)
  var input: File = _

  @Argument(fullName="vep", shortName="vep", doc="path to Variant Effect Predictor perl executable (default: )", required=false)
  var vep: File = _

  @Argument(fullName="python", shortName="python", doc="path to python executable (default: )", required=false)
  var python: File = _

  @Argument(fullName="records_per_job", shortName="records_per_job", doc="number of vcf records to process in each job (default: 1000)", required=false)
  var records_per_job: Int = 1000

  @Argument(fullName="keep_temp_files", shortName="keep_temp_files", doc="keep temporary files", required=false)
  var keep_temp_files: Boolean = false

  override def commandLine = super.commandLine +
                             required (python) +
                             required ("--input", input)
                             optional ("--vep", vep)
  optional ("--vep", vep)


}

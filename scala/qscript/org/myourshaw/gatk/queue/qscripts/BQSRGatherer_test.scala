import org.broadinstitute.sting.queue.function.scattergather

import org.broadinstitute.sting.commandline.Gatherer
import org.broadinstitute.sting.queue.function.InProcessFunction
import collection.JavaConversions._
import java.io.File
import org.broadinstitute.sting.commandline.{Input, Output}
import org.broadinstitute.sting.queue.function.QFunction
import org.broadinstitute.sting.queue.QException
import org.broadinstitute.sting.utils.io.IOUtils
import collection.JavaConversions._

@Output(doc="The original output of the scattered function")
var originalOutput = List("/scratch1/tmp/myourshaw/mmjj_20130514/bams/sample_bams/CIPO10A.sample.recal_data.grp")

@Input(doc="Parts to gather back into the original output")
var gatherParts = List("/scratch1/tmp/myourshaw/mmjj_20130514/sg_dir/.qlog/scratch1/tmp/myourshaw/mmjj_20130514/bams/sample_bams/CIPO10A.sample.recal_data.grp.baseRecalibrator-sg/temp_01_of_86/CIPO10A.sample.recal_data.grp", "/scratch1/tmp/myourshaw/mmjj_20130514/sg_dir/.qlog/scratch1/tmp/myourshaw/mmjj_20130514/bams/sample_bams/CIPO10A.sample.recal_data.grp.baseRecalibrator-sg/temp_02_of_86/CIPO10A.sample.recal_data.grp")


org.broadinstitute.sting.queue.function.scattergather.GathererFunction.run(gatherParts, originalOutput)
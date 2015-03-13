#!/bin/bash

#SGE job array runner
#input: $1 is a job array file with the commands to be run, one per line
#qsub -t or drmaa runBulkJobs calls this script once for each job
#with ${SGE_TASK_ID} = 1-based line number in the job array file

eval "`head -n ${SGE_TASK_ID} ${1} | tail -n 1`";

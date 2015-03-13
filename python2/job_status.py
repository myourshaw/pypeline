#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import time
import datetime
from math import modf
from os.path import isfile

#-a /scratch1/tmp/myourshaw/tmp/accounting -j 1000
#-a /scratch1/tmp/myourshaw/tmp/accounting -j 1028822
#-d 2012-01-21

def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'get SGE accounting info for one or more jobs, optionally limited by user and/or dates',
        epilog = 'pypeline.job_status version 1.0β1 ©2011-2012 Michael Yourshaw all rights reserved')
    parser.add_argument('--accounting', '-a', '-i',
        default='/opt/gridengine/default/common/accounting',
        help='accounting file, usually on rocks54 (default: /opt/gridengine/default/common/accounting)',)
    parser.add_argument('--jobs', '-j', nargs='+', type=int,
        help='list of job ids, or all jobs if omitted',)
    parser.add_argument('--user', '-u',
        help='limit output to jobs with this job owner user ID',)
    parser.add_argument('--dates', '-d', nargs='+',
        help='limit output to date or range of dates (YYYY-mm-dd [YYYY-mm-dd])',)
    parser.add_argument('--warn', '-w', action='store_true', default=False,
        help='warn of corrupt records in accounting file (default: skip without warning)',)
    
    args = parser.parse_args()
    
    first_date, last_date = datetime.date.timetuple(datetime.date.min), datetime.date.timetuple(datetime.date.max)
    if args.dates:
        first_date = time.strptime(args.dates[0], '%Y-%m-%d')
        last_date = time.strptime(args.dates[1]+' 23:59:59', '%Y-%m-%d %H:%M:%S') if len(args.dates)>1 else time.strptime(args.dates[0]+' 23:59:59', '%Y-%m-%d %H:%M:%S')
    first_date = time.mktime(first_date)
    last_date = time.mktime(last_date)
    
    column_doc =\
"""##qname	Name of the cluster queue in which the job has run.
##hostname	Name of the execution host.
##group	The effective group id of the job owner when executing the job.
##owner	Owner of the Sun Grid Engine job.
##job_name	Job name.
##job_number	Job identifier - job number.
##account	An account string as specified by the qsub(1) or qalter(1) -A option.
##priority	Priority value assigned to the job corresponding to the priority parameter in the queue configuration (see queue_conf(5)).
##submission_time	Submission time (GMT unix time stamp).
##start_time	Start time (GMT unix time stamp).
##end_time	End time (GMT unix time stamp).
##failed	Indicates the problem which occurred in case a job could not be started on the execution host (e.g. because the owner of the job did not have a valid account on that machine). If Sun Grid Engine tries to start a job multiple times, this may lead to multiple  entries  in  the  accounting file corresponding to the same job ID.
##exit_status	Exit  status  of  the job script (or Sun Grid Engine specific status in case of certain error conditions).  The exit status  is  determined  by following  the  normal  shell  conventions.   If the command terminates normally the value of the command is its exit status.  However, in  the ase  that  the  command exits abnormally, a value of 0200 (octal), 128 (decimal) is added to the value of the command  to  make  up  the  exit tatus. For  example,  If a job dies through signal 9 (SIGKILL) then the exit status becomes 128 + 9 = 137.
##ru_wallclock	Difference between end_time and start_time (see above).
##ru_utime	This  is  the total amount of time spent executing in user mode, expressed in a timeval structure (seconds plus microseconds).
##ru_stime	This is the total amount of time spent executing in kernel mode, expressed in a timeval structure (seconds plus microseconds).
##ru_maxrss	This  is  the maximum resident set size used (in kilobytes). For RUSAGE_CHILDREN, this is the resident set size  of  the  largest child, not the maximum resident set size of the process tree.
##ru_ixrss	integral shared text size. This field is currently unused on Linux.
##ru_ismrss	integral shared memory size. This field is currently unused on Linux(?).
##ru_idrss	integral unshared data size. This field is currently unused on Linux.
##ru_isrss	integral unshared stack size. This field is currently unused on Linux.
##ru_minflt	The  number  of  page  faults serviced without any I/O activity; here I/O activity is avoided by “reclaiming” a page  frame  from the list of pages awaiting reallocation.
##ru_majflt	The number of page faults serviced that required I/O activity.
##ru_nswap	This field is currently unused on Linux.
##ru_inblock	The number of times the file system had to perform input.
##ru_oublock	The number of times the file system had to perform output.
##ru_msgsnd	messages sent. This field is currently unused on Linux.
##ru_msgrcv	messages received. This field is currently unused on Linux.
##ru_nsignals	signals received. This field is currently unused on Linux.
##ru_nvcsw	The  number  of times a context switch resulted due to a process voluntarily giving up the processor before its  time  slice  was completed (usually to await availability of a resource).
##ru_nivcsw	The  number  of  times a context switch resulted due to a higher priority  process  becoming  runnable  or  because  the  current process exceeded its time slice.
##project	The project which was assigned to the job.
##department	The department which was assigned to the job.
##granted_pe	The parallel environment which was selected for that job.
##slots	The number of slots which were dispatched to the job by the scheduler.
##task_number	Array job task index number.
##cpu	The cpu time usage in seconds.
##mem	The integral memory usage in Gbytes cpu seconds.
##io	The amount of data transferred in input/output operations.
##category	A string specifying the job category.
##iow	The io wait time in seconds.
##pe_taskid	If  this  identifier is set the task was part of a parallel job and was passed to Sun Grid Engine via the qrsh -inherit interface.
##maxvmem	The maximum vmem size in bytes.
##arid	Advance reservation identifier. If the job used resources of an advance reservation  then  this  field  contains  a positive integer identifier otherwise the value is "0" .
##ar_submission_time	If the job used resources of an advance  reservation  then  this  field contains  the  submission  time  (GMT  unix  time stamp) of the advance reservation, otherwise the value is "0" ."""
    columns = ['record_no','qname','hostname','group','owner','job_name','job_number','account','priority','submission_time','start_time','end_time','failed','exit_status','ru_wallclock','ru_utime','ru_stime','ru_maxrss','ru_ixrss','ru_ismrss','ru_idrss','ru_isrss','ru_minflt','ru_majflt','ru_nswap','ru_inblock','ru_oublock','ru_msgsnd','ru_msgrcv','ru_nsignals','ru_nvcsw','ru_nivcsw','project','department','granted_pe','slots','task_number','cpu','mem','io','category','iow','pe_taskid','maxvmem','arid','ar_submission_time']

    all_jobs = not bool(args.jobs)
    if args.jobs:
        args.jobs = sorted(set(args.jobs))
    line_count = 0
    
    try:
        with open(args.accounting,'r') as acct:
            print column_doc
            print '#'+'\t'.join(columns)
            for line in acct:
                line_count+=1
                if line.startswith('#') or line.strip() == '':
                    continue
                line = line.rstrip('\n')
                fields = line.split(':')
                fields.insert(0,str(line_count))
                if len(fields) != len(columns):
                    if args.warn:
                        sys.stderr.write('WARNING CORRUPT LINE {}:\n{}\n'.format(line_count,line))
                    continue
                record = {columns[i]: fields[i] for i in range(len(fields))}
                if ((all_jobs or (args.jobs and int(record['job_number']) in args.jobs)) \
                    and (not args.user or record['owner']==args.user) \
                    and ((float(record['submission_time'])>=first_date and float(record['submission_time'])<=last_date) or(float(record['start_time'])>=first_date and float(record['start_time'])<=last_date) or (float(record['end_time'])>=first_date and float(record['end_time'])<=last_date))):
                    record['submission_time'] = time.strftime('%Y-%m-%d %H:%M', time.localtime(int(record['submission_time'])))
                    record['start_time'] = time.strftime('%Y-%m-%d %H:%M', time.localtime(int(record['start_time'])))
                    record['end_time'] = time.strftime('%Y-%m-%d %H:%M', time.localtime(int(record['end_time'])))
                    record['ru_wallclock'] = str(datetime.timedelta(seconds=int(record['ru_wallclock'])))
                    record['ru_utime'] = str(datetime.timedelta(seconds=float(record['ru_utime'])))
                    record['ru_stime'] = str(datetime.timedelta(seconds=float(record['ru_stime'])))
                    record['cpu'] = str(datetime.timedelta(seconds=float(record['cpu'])))
                    line_out = '\t'.join([record[columns[i]] for i in range(len(columns))])
                    print line_out
    except IOError:
        return 'Could not open {}.\nThis is usually caused by running this application on a cluster node.\nTry running it on rocks.'.format(args.accounting)
    else:
        return 'done'
        
if __name__ == "__main__": sys.exit(main())

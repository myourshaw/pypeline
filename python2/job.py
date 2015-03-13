#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import stat
import tempfile
import shutil
import re
import math
try:
 import drmaa
except ImportError, e:
 pass # module doesn't exist, deal with it. pip install drmaa
import my

class JobError(Exception): pass
class JobTerminatedAbnormallyException(Exception): pass
class SessionExistsException(Exception): pass

def makedir(dir):
	try:
	  os.makedirs(dir)
	except OSError:
		if os.path.isdir(dir):
			pass
		else:
			raise

jobNameBadCharsRegex = re.compile(r"[\s/:@\\\*\?]")

jobStatusHeader = 'jobId\tjobName\twasAborted\thasExited\texitStatus\thasSignal\tterminatedSignal\thasCoreDump\tresourceUsage'

def format_status_msg(jobname, status):
	return '{jobId}\t{0}\t{wasAborted}\t{hasExited}\t{exitStatus}\t{hasSignal}\t{terminatedSignal}\t{hasCoreDump}\t{resourceUsage}'\
	.format(jobname, **status._asdict())

class Job(object):
    def __init__(self,
        remoteCommand = None,
        args = None,
        jobName = None,
        memory = None,
        processors = 1,
        queue = None,
        hold_jid = None,
        queuingModeBunchJobsOnNodes = True,
        email = None,
        nativeSpecification = None,
        inputPath = None,
        outputPath = None,
        errorPath = None,
        workingDirectory = None,
        blockEmail = None,
        jobEnvironment = None,
        startTime = None,
        deadlineTime = None,
        softRunDurationLimit = None,
        hardRunDurationLimit = None,
        softWallclockTimeLimit = None,
        hardWallclockTimeLimit = None,
        jobCategory = None,
        jobSubmissionState = None,
        joinFiles = None,
        transferFiles = None,
        exclusive = True,
        ):
        self.remoteCommand = remoteCommand
        self.args = args
        self.jobName = jobName
        self.memory = memory
        self.processors = processors
        self.queue = queue
        self.hold_jid = hold_jid
        self.queuingModeBunchJobsOnNodes = queuingModeBunchJobsOnNodes
        self.nativeSpecification = nativeSpecification
        self.email = email
        self.inputPath = inputPath
        self.outputPath = outputPath
        self.errorPath = errorPath
        self.workingDirectory = workingDirectory
        self.blockEmail = blockEmail
        self.jobEnvironment = jobEnvironment
        self.startTime = startTime
        self.deadlineTime = deadlineTime
        self.softRunDurationLimit = softRunDurationLimit
        self.hardRunDurationLimit = hardRunDurationLimit
        self.softWallclockTimeLimit = softWallclockTimeLimit
        self.hardWallclockTimeLimit = hardWallclockTimeLimit
        self.jobCategory = jobCategory
        self.jobSubmissionState = jobSubmissionState
        self.joinFiles = joinFiles
        self.transferFiles = transferFiles
        self.exclusive = exclusive
        
        self.session = None
        self.jobTemplate = None
        self.jobId = None
        self.jobList = None
        self.completionStatus = None
        self.completionMsg = None

    def setup_job(self):
        if self.session:
            raise SessionExistsException('A drmaa sesion already exists for this job')
        self.session = drmaa.Session()
        self.session.initialize()
        self.jobTemplate = self.session.createJobTemplate()
        if self.remoteCommand: self.jobTemplate.remoteCommand = self.remoteCommand
        if self.args: self.jobTemplate.args = self.args
        self.jobTemplate.jobName = self.jobName if self.jobName else os.path.basename(self.jobTemplate.remoteCommand)
        self.jobTemplate.jobName = jobNameBadCharsRegex.sub('_', self.jobTemplate.jobName)
        __spec = [self.nativeSpecification] if self.nativeSpecification else ['-V']
        #__memqueue = 'all.q@compute-4*'
        __memqueue = None
        if self.memory:
            __memMG = self.memory.upper()[-1]
            __memx = pow(2,20) if __memMG == 'M' else pow(2,30) if __memMG == 'G' else 1
            __mem = math.ceil(__memx * float(self.memory.upper()[:-1]))
            if __mem > 7* pow(2,30): __memqueue = 'all.q@compute-[235]*'
        if __memqueue or self.queue:
            __spec.append('-q {}'.format(self.queue if self.queue else __memqueue))
        #if self.memory: __spec.append('-hard -l mem_free={}'.format(self.memory))
        #if self.memory: __spec.append('-l mem_free={}'.format(self.memory))
        if self.processors:
            if int(self.processors) > 1:
                __spec.append('-hard -pe serial {}'.format(self.processors))
            else:
                __spec.append('-hard -pe {} 1'.format('make' if self.queuingModeBunchJobsOnNodes and not self.memory else 'serial'))
        if self.exclusive:
            __spec.append('-w n -l excl=true')
        if self.hold_jid:
            if isinstance(self.hold_jid, (tuple, list)):
                __spec.append('-hold_jid {}'.format(','.join(self.hold_jid)))
            elif isinstance(self.hold_jid, (basestring, int)):
                __spec.append('-hold_jid {}'.format(self.hold_jid))
            else:
                raise JobError('incorrect type {} for hold_jid {}'.format(type(hold_jid), hold_jid))
        self.jobTemplate.nativeSpecification = ' '.join(__spec)
        if self.email: self.jobTemplate.email = self.email
        if self.inputPath and not self.inputPath.startswith(':'):
            self.inputPath = ':' + self.inputPath
        if self.inputPath: self.jobTemplate.inputPath = self.inputPath
        if self.outputPath and not self.outputPath.startswith(':'):
            self.outputPath = ':' + self.outputPath
        if self.outputPath: self.jobTemplate.outputPath = self.outputPath
        if self.errorPath and not self.errorPath.startswith(':'):
            self.errorPath = ':' + self.errorPath
        if self.errorPath: self.jobTemplate.errorPath = self.errorPath
        if self.workingDirectory: self.jobTemplate.workingDirectory = self.workingDirectory
        if self.blockEmail != None: self.jobTemplate.blockEmail = self.blockEmail
        if self.jobEnvironment: self.jobTemplate.jobEnvironment = self.jobEnvironment
        if self.startTime: self.jobTemplate.startTime = self.startTime
        if self.deadlineTime: self.jobTemplate.deadlineTime = self.deadlineTime
        if self.softRunDurationLimit: self.jobTemplate.softRunDurationLimit = self.softRunDurationLimit
        if self.hardRunDurationLimit: self.jobTemplate.hardRunDurationLimit = self.hardRunDurationLimit
        if self.softWallclockTimeLimit: self.jobTemplate.softWallclockTimeLimit = self.softWallclockTimeLimit
        if self.hardWallclockTimeLimit: self.jobTemplate.hardWallclockTimeLimit = self.hardWallclockTimeLimit
        if self.jobCategory: self.jobTemplate.jobCategory = self.jobCategory
        if self.jobSubmissionState: self.jobTemplate.jobSubmissionState = self.jobSubmissionState
        if self.joinFiles: self.jobTemplate.joinFiles = self.joinFiles
        if self.transferFiles: self.jobTemplate.transferFiles = self.transferFiles

    def executeShellScript(self, script, args, synchronous=True):
        self.remoteCommand = script
        self.args = args
        self.setup_job()
        try:
            self.jobId = self.session.runJob(self.jobTemplate)
            if synchronous:
                self.completionStatus = __status = self.session.wait(self.jobId, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        except Exception as __e:
            raise __e
        else:
            if synchronous:
                self.completionMsg = [self.jobId, jobStatusHeader]
                __failed_jobs = [self.jobId, jobStatusHeader]
                self.completionMsg.append(format_status_msg(self.jobTemplate.jobName, __status))
                if self.completionStatus.wasAborted or not self.completionStatus.hasExited or self.completionStatus.exitStatus:
                    __failed_jobs = [__failed_jobs.append(format_status_msg(self.jobTemplate.jobName, __status))]
                    raise JobTerminatedAbnormallyException('job terminated abnormally\n{}\nStatus is:\n{}\n'.format('\n'.join(__failed_jobs), '\n'.join(self.completionMsg)))
            else:
                self.completionMsg = [self.jobId, 'completion status not available for asynchronous job']
        finally:
            self.session.deleteJobTemplate(self.jobTemplate)
            self.session.exit()

    def executeCommand(self, command, synchronous=True):
        #tmp dir and files for commands and shell script
        __prefix = 'job_{}_tmp'.format(self.jobName if self.jobName else '')
        __dir = self.workingDirectory if self.workingDirectory else os.path.expanduser('~')
        self.tmp_dir = my.unique_dir(__prefix, __dir)
        __script_file, __script_filename = my.unique_file(os.path.join(self.tmp_dir, __prefix+'.sh'))
        __script_file.write("""#!/bin/bash
set -o pipefail;
eval "${1}";
rc=$?
if [[ $rc != 0 ]]; then
  exit 100;
fi
""")
        __script_file.close()
        os.chmod(__script_filename, stat.S_IRWXU)
        self.remoteCommand = __script_filename
        self.args = [command]
        self.setup_job()
        try:
            self.jobId = self.session.runJob(self.jobTemplate)
            if synchronous:
                self.completionStatus = __status = self.session.wait(self.jobId, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        except Exception as __e:
            raise __e
        else:
            if synchronous:
                self.completionMsg = [self.jobId, jobStatusHeader]
                __failed_jobs = [self.jobId, jobStatusHeader]
                self.completionMsg.append(format_status_msg(self.jobTemplate.jobName, __status))
                if self.completionStatus.wasAborted or not self.completionStatus.hasExited or self.completionStatus.exitStatus:
                    __failed_jobs = [__failed_jobs.append(format_status_msg(self.jobTemplate.jobName, __status))]
                    raise JobTerminatedAbnormallyException('job terminated abnormally\n{}\nStatus is:\n{}\n'.format('\n'.join(__failed_jobs), '\n'.join(self.completionMsg)))
            else:
                self.completionMsg = [self.jobId, 'completion status not available for asynchronous job']
        finally:
            self.session.deleteJobTemplate(self.jobTemplate)
            self.session.exit()

    def executeCommands(self, commands, synchronous=True):
        if isinstance(commands, (list, tuple)):
            array_job = True
        elif isinstance(commands, (basestring)):
            array_job = False
        else:
            raise JobError('bad type {} for command(s) {}'.format(type(commands), commands))
        #tmp dir and files for commands and shell script
        __prefix = 'job_{}_tmp'.format(self.jobName if self.jobName else '')
        __dir = self.workingDirectory if self.workingDirectory else os.path.expanduser('~')
        self.tmp_dir = my.unique_dir(__prefix, __dir)
        if array_job:
            __command_file, __command_filename = my.unique_file(os.path.join(self.tmp_dir, __prefix+'.commands'))
            __command_file.writelines("%s\n" % item for item in commands)
            __command_file.close()
            __script_file, __script_filename = my.unique_file(os.path.join(self.tmp_dir, __prefix+'.sh'))
            __script_file.write("""#!/bin/bash
set -o pipefail;
eval "`head -n ${SGE_TASK_ID} ${1} | tail -n 1`";
rc=$?
if [[ $rc != 0 ]]; then
  exit 100;
fi
""")
            self.args = [__command_filename]
        else:
            __script_file, __script_filename = my.unique_file(os.path.join(self.tmp_dir, __prefix+'.sh'))
            __script_file.write("""#!/bin/bash
set -o pipefail;
eval "${1}";
rc=$?;
if [[ $rc != 0 ]]; then
  exit 100;
fi
""")
            self.args = [commands]
        __script_file.close()
        os.chmod(__script_filename, stat.S_IRWXU)
        self.remoteCommand = __script_filename
        self.setup_job()
        try:
            if array_job:
                self.jobList = self.session.runBulkJobs(self.jobTemplate, 1, len(commands), 1)
                self.jobId = self.jobList[0].split('.')[0]
                if synchronous:
                    self.session.synchronize(self.jobList, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
            else:
                self.jobId = self.session.runJob(self.jobTemplate)
                if synchronous:
                    self.completionStatus = __status = self.session.wait(self.jobId, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        except Exception as __e:
            raise __e
        else:
            if array_job:
                if synchronous:
                    self.completionStatus = []
                    self.completionMsg = [self.jobId, jobStatusHeader]
                    for __j in self.jobList:
                        self.completionStatus.append(self.session.wait(__j, drmaa.Session.TIMEOUT_WAIT_FOREVER))
                    __jobs_failed = False
                    __failed_jobs = [jobStatusHeader]
                    for __status in self.completionStatus:
                        self.completionMsg.append(format_status_msg(self.jobTemplate.jobName, __status))
                        if __status.wasAborted or not __status.hasExited or __status.exitStatus:
                            __jobs_failed = True
                            __failed_jobs.append(format_status_msg(self.jobTemplate.jobName, __status))
                    if __jobs_failed:
                        raise JobTerminatedAbnormallyException('One or more commands terminated abnormally\n{}\nStatus of all jobs is:\n{}\n'.format('\n'.join(__failed_jobs), '\n'.join(self.completionMsg)))
                    shutil.rmtree(self.tmp_dir, ignore_errors=True)
                else:
                    self.completionMsg = [self.jobId, 'completion status not available for asynchronous jobs']
            else:
                if synchronous:
                    self.completionMsg = [self.jobId, jobStatusHeader]
                    __failed_jobs = [jobStatusHeader]
                    self.completionMsg.append(format_status_msg(self.jobTemplate.jobName, __status))
                    if self.completionStatus.wasAborted or not self.completionStatus.hasExited or self.completionStatus.exitStatus:
                        __failed_jobs = [__failed_jobs.append(format_status_msg(self.jobTemplate.jobName, __status))]
                        raise JobTerminatedAbnormallyException('job terminated abnormally\n{}\nStatus is:\n{}\n'.format('\n'.join(__failed_jobs), '\n'.join(self.completionMsg)))
                    shutil.rmtree(self.tmp_dir, ignore_errors=True)
                else:
                    self.completionMsg = [self.jobId, 'completion status not available for asynchronous job']
        finally:
            self.session.deleteJobTemplate(self.jobTemplate)
            self.session.exit()

    def suspend(self):
        self.session.control(self.jobId, drmaa.JobControlAction.SUSPEND)
    
    def resume(self):
        self.session.control(self.jobId, drmaa.JobControlAction.RESUME)
    
    def hold(self):
        self.session.control(self.jobId, drmaa.JobControlAction.HOLD)
    
    def release(self):
        self.session.control(self.jobId, drmaa.JobControlAction.RELEASE)
    
    def terminate(self):
        self.session.control(self.jobId, drmaa.JobControlAction.TERMINATE)

    @property
    def jobStatus(self):
        decodestatus = {
            drmaa.JobState.UNDETERMINED: 'UNDETERMINED',
            drmaa.JobState.QUEUED_ACTIVE: 'QUEUED_ACTIVE',
            drmaa.JobState.SYSTEM_ON_HOLD: 'SYSTEM_ON_HOLD',
            drmaa.JobState.USER_ON_HOLD: 'USER_ON_HOLD',
            drmaa.JobState.USER_SYSTEM_ON_HOLD: 'USER_SYSTEM_ON_HOLD',
            drmaa.JobState.RUNNING: 'RUNNING',
            drmaa.JobState.SYSTEM_SUSPENDED: 'SYSTEM_SUSPENDED',
            drmaa.JobState.USER_SUSPENDED: 'USER_SUSPENDED',
            drmaa.JobState.DONE: 'DONE',
            drmaa.JobState.FAILED: 'FAILED',
            }
        self.session.initialize(self.session.contact)
        return (decodestatus(self.session.jobStatus(self.jobId)))

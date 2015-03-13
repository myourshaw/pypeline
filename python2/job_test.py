#!/usr/bin/env python

import sys
import os
import stat
import tempfile
import job

synchronous = True
commands = 2
job = job.Job()
if commands:
	cmds = []
	for i in range(1,commands+1,2):
		cmds.append('{} {} {}'.format(os.path.join(os.getcwd(),'sleeper.sh'),os.environ['LOGNAME'],str(i)))
	job.executeCommands(cmds, synchronous)
else:
	job.executeShellScript(os.path.join(os.getcwd(),'sleeper.sh'),[os.environ['LOGNAME'],'30'], synchronous)
print '\n'.join(job.completionMsg)
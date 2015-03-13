#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from ConfigParser import SafeConfigParser #configparser in python 3
import shutil
import glob
import my
import job

dir = '/home/myourshaw/foo'
my.makedir(dir)
cmds = []
for i in range(2):
	cmds.append('sleep 1')
j = my.run_job(cmds, jobName='sleep', job_dir=dir, qout_dir=dir, email=None, synchronous=True, processors=None, memory=None, hold_jid=None)
id = j.jobId


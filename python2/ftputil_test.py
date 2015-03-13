#!/usr/bin/env python

import os
from ftputil import FTPHost
host = 'ftp.omim.org'
remote_directory = '/OMIM/'
local_directory = '/home/myourshaw/tmp'
with FTPHost(host, 'anonymous', '') as h:
    h.chdir(remote_directory)
    names = h.listdir(h.curdir)
    for name in names:
        if h.path.isfile(name):
            h.download(name, os.path.join(local_directory, name))
            print(name)

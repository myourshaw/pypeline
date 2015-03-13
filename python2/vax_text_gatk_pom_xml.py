#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import my


tree = ET.ElementTree(file='/share/apps/myourshaw/gatk_2014-05-02/pom.xml')
root = tree.getroot()
print root.tag, root.attrib
#http://stackoverflow.com/questions/1953761/accessing-xmlns-attribute-with-python-elementree
if root.tag[0] == '{':
    uri, ignore, tag = root.tag[1:].partition('}')
else:
    uri = None
    tag = root.tag
for node in root:
    print node.tag, node.attrib
ns_parent = 'parent' if uri == None else '{{{}}}parent'.format(uri)
ns_version = 'parent/version' if uri == None else '{{{0}}}parent/{{{0}}}version'.format(uri)
parent_node = root.find(ns_parent)
version_node = root.find(ns_version)
version = version_node.text
print version



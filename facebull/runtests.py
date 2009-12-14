#!/bin/python

import sys
import os
import dircache
import time
import filecmp

testdir = "./test"
in_ext = '.in'
out_ext = '.out'
run_ext = '.run'
filelist = dircache.listdir(testdir)

#get ride of duplicates
testnames = list(set(map(lambda x: x.split('.')[0], filelist)))
testnames.sort()

for test in testnames:
    t0 = time.time()
    infile = os.path.join(testdir, test+in_ext)
    runfile = os.path.join(testdir,test+run_ext)
    os.system("./facebull " + infile + " >" + runfile)
    t1 = time.time()

    outfile = os.path.join(testdir, test+out_ext)
    f1 = open(runfile).read()
    f2 = open(outfile).read()
    if f1 == f2:
        print "%s -> PASS (%0.4fs)" % (test, t1-t0)
    else:
        print "\t%s -> FAIL (%0.4fs)" % (test, t1-t0)
    


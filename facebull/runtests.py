#!/bin/python

import sys
import os
import dircache
import time
import subprocess
import filecmp
from threading import Thread
testpid = -1

class TestRunner(Thread):
    def __init__(self, tests=[], binaries=[]):
        Thread.__init__(self)
        self._tests = tests
        self._binaries = binaries
        return

    def run(self):
        global testpid
        testpid = os.getpid()
        try:
            for test in self._tests:
                for bin in self._binaries:
                    t0 = time.time()
                    infile = os.path.join(testdir, test+in_ext)
                    runfile = os.path.join(testdir,test+run_ext)    
                    
                    f1 = open(runfile, "w+")
                    subprocess.check_call([bin, infile], stdout=f1)
                    #subprocess.call([bin, infile], stdout=f1)                    
                    t1 = time.time()
                    outfile = os.path.join(testdir, test+out_ext)
                    f1.close()
                    f1 = open(runfile).read()
                    f2 = open(outfile).read()
                    if f1 == f2:
                        print "\t%s -> PASS (%0.5f s)" % (bin[2:], t1-t0)
                    else:
                        print "\t%s -> FAIL (%0.5f s)" % (bin[2:], t1-t0)
        except KeyboardInterrupt:
            #    print "thread.run(): Recevied Ctrl-C"
            return


def TestSuper(tests, algs, timeout=60):
    try:
        for test in tests:
            print "\nRUNNING TEST %s" % test
            for alg in algs:
                tester = TestRunner([test], [alg])
                tester.start()#run()
                tester.join(timeout)
                if tester.isAlive():
                    print "test timeout"
                    os.popen("kill -9 "+str(testpid))
                    return

    except KeyboardInterrupt:
        print ""#Received Ctrl-C in Main Thread."



testdir = "./test"
in_ext = '.in'
out_ext = '.out'
run_ext = '.run'
def main():
    filelist = dircache.listdir(testdir)
    timeout = 60

    #get ride of duplicates
    testnames = list(set(map(lambda x: x.split('.')[0], filter(lambda x: x.endswith('.in'), filelist))))

    testnames.sort()
    testnames.reverse()
    algorithms = ['./pathwalk','./edgewalk','./megbb']
    tester = None

    if len(sys.argv) == 2:
        if sys.argv[1] == 'pathwalk':
            TestSuper(testnames, [algorithms[0]], timeout)
        elif sys.argv[1] == 'edgewalk':
            TestSuper(testnames, [algorithms[1]], timeout)        
        elif sys.argv[1] == 'megbb':
            TestSuper(testnames, [algorithms[2]], timeout)
        else:
            tester = TestSuper([sys.argv[1]], algorithms, timeout)
    else:
        tester = TestSuper(testnames, algorithms, timeout)

    sys.exit(0)
        

if __name__ == '__main__':
    main()


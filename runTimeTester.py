#This script outputs a csv of run time results for different grid sizes n.

from subprocess import call
import os

for n in range(32, 71):
    print "n is %d" % (n)
    os.system("./gs_solver --grid-elems-r %d --grid-elems-z %d | tail -n 1 | awk '{ print \"%d,\", $4, $8}' >> runTimeResults.txt" % (n,n,n))

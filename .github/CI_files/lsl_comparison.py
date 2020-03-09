#!/usr/bin/env python
# A script used within CI to check if the results of lsl have changed.

import sys
file1 = open(sys.argv[1],'r')
file2 = open(sys.argv[2],'r')
eps = 0.00000000001

for i in file1:
    if(i[0] == "#"):
        pass
    else:
        f1_results = i.split()

for i in file2:
    if(i[0] == "#"):
        pass
    else:
        f2_results = i.split()

for i in range(len(f1_results)):
    if (i == 0):
        pass
    else:
        assert(abs(float(f1_results[i]) - float(f2_results[i])) < eps)



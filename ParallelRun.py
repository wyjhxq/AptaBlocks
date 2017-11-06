import random
import math
import subprocess
import os
import numpy

#==================================================================  files
Lambda1 = arange(0.5, 10.5, 0.5) # energy
Lambda2 = arange(0.5, 1.7, 0.2)

print Lambda1
print Lambda2

IterNum = 3

commandlist = []
for iter in range(IterNum):
    for lbda in Lambda1:
        for lbda2 in Lambda2:
#            outfilename = "1102out_L1-"+str(lbda)+"&L2-"+str(lbda2)+"_"+str(iter)+"_Shasha_21mer_spacer.txt"
#commandlist.append([ "./AptaBlockAptamer_siRNA", str(lbda), str(lbda2), str(outfilename)])


#processes = set()
#max_processes = 30
#for i in range(len(commandlist)):
#    print(commandlist[i])
#    processes.add(subprocess.Popen(commandlist[i]))
#    if len(processes) >= max_processes:
#        os.wait()
#        processes.difference_update([p for p in processes if p.poll() is not None])










































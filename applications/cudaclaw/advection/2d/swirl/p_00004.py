# comment = "Torus example : eff. resolution = 2048 x 2048"
import sys
import os
import subprocess
import random

np = 4
arg_list = ["mpirun","-n",str(np),"swirl"]
jobid = random.randint(1000,9999)
outfile = "swirl_0000{:d}.o{:d}".format(np,jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process {:d} with jobid {:d} on {:d} processor(s).".format(po.pid,jobid,np))
# po.wait()

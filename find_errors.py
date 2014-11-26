from os import listdir, system
from os.path import isfile, join
#sorts through all of the stupid .o######## files the HPC outputs
#deletes them if they are unimportant, reads errors otherwise

mypath="./" #directory path, i'm already here
prefix = 'IMSRG'  #five letter prefix for ".o" output files

files = [ f for f in listdir(mypath) if isfile(join(mypath,f))]
outfiles=[]

for f in files:
    if f[:5] == prefix:
        outfiles.append(f) 

for fname in outfiles: 
   # delete = True
    fx = open(fname,'r') 

    for line in fx:
        if "FAIL" in line:
            print 'secondary convergence fail: ', fname[:,-10]
        if "walltime" in line:
            print 'walltime exceeded: ', fname[:-10]
            delete=False
            break
        if "mem" in line:
            print 'memory exceeded: ', fname[:-10]
            delete=False
            break
        if "error" in line: 
            print 'some stupid error: ', fname[:-10]
            delete=False
            break
        if "fault" in line:
            print 'segmentation fault: ', fname[:-10]
            delete=False
            break
        if "Job" in line:
            break
   #if delete:
     #   system('rm '+fname) 
        


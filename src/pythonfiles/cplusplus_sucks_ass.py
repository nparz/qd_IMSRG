# because I don't want to deal with C++ IO stuff. 
# I'm a wimp
import sys

s = sys.argv[1].strip()
eoff = float(sys.argv[2].strip())
fl = open('evals.dat','r')

outstr = s + '  ' 
for line in fl:
    outx = str(eoff + float(line.strip()))
    outstr = outstr + outx + '  '

outstr = outstr+'\n'
fl.close()

fx = open('flowing_CI.dat','a') 

fx.write(outstr)
fx.close()
    

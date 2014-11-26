from math import exp,log
import os

nstr = raw_input('enter amount of electrons, seperated by commas: ')
nstr=nstr.strip()
nlist = nstr.split(',')

hwstr = raw_input('enter hw values, seperated by commas: ')
hwstr=hwstr.strip()
hwlist = hwstr.split(',') 

Rstr= raw_input('enter number of shells, seperated by commas: ')
Rstr=Rstr.strip()
Rlist=Rstr.split(',')

ty = raw_input( 'enter M for magnus, T for traditional or Q for matrix elements: ')  


wtimes = [ 50 , 400, 1000, 4000, 8000, 12000, 14350,25000,35000] 

if ty == 'M':
    
    fq = open('run_mag.bat','w')

    fq.write( '#!/bin/bash \n\n')
    
    for n in nlist:
        for hw in hwlist:
            for R in Rlist: 
                nx = float(n)
                hwx = float(hw) 
                Rx= float(R) 
                try: 
                    
                    if '.' not in hw:
                        hw=hw+'.'
    
                    carg = n+' '+hw+'d0'+' '+R+'\n\n'

                except IndexError:
                    print "Write Failed: incorrect number of arguments"
                    raise SystemExit
                
       
                if hw[0] == '.':
                    hw = '0'+hw
                while len(hw) < 4:
                    hw = hw +'0'
                
                fq.write( 'qsub pbs_'+n+'_'+hw+'_'+R+'_magnus\n') 
                
                fl = open('pbs_'+n+'_'+hw+'_'+R+'_magnus','w')
                
                jobname = 'IMSRG_'+n+'_'+hw+'_'+R+'_magnus'
                if R == '3' or R=='4':
                    mem = '10mb'
                elif R=='5' or R=='6': 
                    mem = '100mb' 
                elif R=='7' or R=='8':
                    mem = '500mb'
                elif R == '9':
                    mem = '1gb'
                else:
                    mem = '2gb'

 #               fac = log(nx) - log(.73575)
#                wall = int(exp( -11.0 + 4.6*log(Rx*Rx+Rx) )*fac/hwx/5.0) 
                    
            
                wall = wtimes[int(R)-3]
                if R == 7:   # density calculated here, so it takes a lot longer
                    wall = wall * 2.0 
                    
                hours = wall/3600 
                wall = wall - hours * 3600
                minutes = wall/60
                wall = wall - minutes * 60
                seconds = wall
                
                if seconds == 0 or seconds == 1: 
                    seconds = 2

                hrs = str(hours)
                if len(hrs) == 1:
                    hrs = '0'+hrs

                mns = str(minutes)
                if len(mns) == 1:
                    mns = '0'+mns
                    
                secs = str(seconds)
                if len(secs) == 1:
                    secs = '0'+secs
                    
                walltime = hrs+':'+mns+':'+secs
                
                nodes = '1'
                cores = '8'
                fl.write( '#!/bin/sh \n\n' ) 
                fl.write( '#PBS -l walltime='+walltime+'\n')
                fl.write( '#PBS -l nodes=1:ppn='+cores+'\n')
                fl.write( '#PBS -l mem='+mem+'\n')
                fl.write( '#PBS -j oe\n')
                fl.write( '#PBS -N '+jobname+'\n')
              # fl.write( '#PBS -M parzuchowski@frib.msu.edu\n')
               # fl.write( '#PBS -m a \n\n') 


                fl.write( 'cd $HOME/qd_imsrg/src/ \n\n')
                fl.write( 'export OMP_NUM_THREADS='+cores+'\n\n' ) 
       
                fl.write('./run_magnus '+carg) 
                fl.write('qstat -f ${PBS_JOBID}')
                fl.close()
	
                os.system("chmod 0755 run_mag.bat" ) 

elif ty == 'T':
    
    gs = raw_input('Is this a ground state calculation (y/n)? ')
    if (gs.lower() == 'y'):
        batchfile='run_ground.bat'
        addlol = ''
    else:
        ml,ms,cut = raw_input('Enter ML , MS ,SHELL CUT: ').split(",")
        batchfile = 'run_excited_'+ml+'_'+ms+'.bat'
        addlol = '_'+ml+'_'+ms
    fq = open(batchfile,'w')
    fq.write( '#!/bin/bash \n\n')
    
    #fk = open('SRC/cases_trad.dat','w') 
    #fk.write('0\n')
    for n in nlist:
        for hw in hwlist:
            for R in Rlist:
                nx = float(n)
                hwx = float(hw) 
                Rx= float(R) 
                try: 
                    
                    if '.' not in hw:
                        hw=hw+'.'
    
                    
                    carg = n+' '+hw+'d0'+' '+R+'\n\n'
                        
                    

                except IndexError:
                    print "Write Failed: incorrect number of arguments"
                    raise SystemExit
                
                
                if hw[0] == '.':
                    hw = '0'+hw
                while len(hw) < 4:
                    hw = hw +'0'

                fq.write( 'qsub pbs_'+n+'_'+hw+'_'+R+addlol+'\n') 
                
                fl = open('pbs_'+n+'_'+hw+'_'+R+addlol,'w')
                
                jobname = 'IMSRG_'+n+'_'+hw+'_'+R
                if R == '3' or R=='4':
                    mem = '10mb'
                elif R=='5' or R=='6': 
                    mem = '100mb' 
                elif R=='7' or R=='8':
                    mem = '500mb'
                else:
                    mem = '1gb' 
               
                if R=='4':
                    mem = '200mb'
                if R=='5':
                    mem = '3gb'
                #if R=='6':
                 #   mem = '10gb'
                #fac = log(nx) - log(.73575)
                #wall = 12*int(2 * exp( -11.0 + 4.3*log(Rx*Rx+Rx) ) * fac / hwx) 
                wall = wtimes[int(R)-3]
                hours = wall/3600 
                wall = wall - hours * 3600
                minutes = wall/60
                wall = wall - minutes * 60
                seconds = wall
                
                if seconds == 0 or seconds == 1: 
                    seconds = 2
                    
                hrs = str(hours)
                if len(hrs) == 1:
                    hrs = '0'+hrs

                mns = str(minutes)
                if len(mns) == 1:
                    mns = '0'+mns
                    
                secs = str(seconds)
                if len(secs) == 1:
                    secs = '0'+secs
                    
                walltime = hrs+':'+mns+':'+secs
                
                nodes = '1'
                cores = '8'
                fl.write( '#!/bin/sh \n\n' ) 
                fl.write( '#PBS -l walltime='+walltime+'\n')
                fl.write( '#PBS -l nodes=1:ppn='+cores+'\n')
                fl.write( '#PBS -l mem='+mem+'\n')
                fl.write( '#PBS -j oe\n')
                fl.write( '#PBS -N '+jobname+'\n')
               # fl.write( '#PBS -M parzuchowski@frib.msu.edu\n')
               # fl.write( '#PBS -m a \n\n') 
            

                fl.write( 'cd $HOME/qd_imsrg/src/ \n\n')
                fl.write( 'export OMP_NUM_THREADS='+cores+'\n\n' ) 
                
                if (gs.lower() == 'y'):
                    
                    fl.write('./gs_decouple ' + carg ) 
                else:
                    carg = carg[:-4]+' '+ml+' '+ms+' '+cut+'\n\n'
                    fl.write('./ex_decouple ' + carg ) 
                
                    
                fl.write('qstat -f ${PBS_JOBID}')
                fl.close()

                os.system("chmod 0755 "+batchfile ) 

elif ty == 'Q':
    
    fq = open('run_TBME.bat','w')

    fq.write( '#!/bin/bash \n\n')

    #fk = open('SRC/cases_TBME.dat','w') 
    #fk.write('0\n')
    for hw in hwlist:
        for R in Rlist: 
               
            hwx = float(hw) 
            Rx= float(R) 
            try:               

                if '.' not in hw:
                    hw=hw+'.'
           
                carg = hw+'d0'+' '+R+'\n\n'

            except IndexError:
                print "Write Failed: incorrect number of arguments"
                raise SystemExit
                
                
            if hw[0] == '.':
                hw = '0'+hw
            while len(hw) < 4:
                hw = hw +'0'

            fq.write( 'qsub pbs_'+hw+'_'+R+'_TBME\n') 
                
            fl = open('pbs_'+hw+'_'+R+'_TBME','w')
                
            jobname = 'IMSRG_'+hw+'_'+R+'_TBME'
        
            memory = int(((Rx*Rx+Rx)*(Rx*Rx+Rx-1))**2*4.25)
                
            if memory/1000000000 != 0:
                mem = '800mb'
                mem = str(memory/1000000000 + 1)+'gb'
            elif  memory/1000000 != 0:
                mem = str(memory/1000000 + 1)+'mb'
            elif  memory/1000 != 0:
                mem = str(memory/1000+1)+ 'kb' 
            else:
                mem = '1kb' 
              
            wall = max(int(exp( -21.0 + 6.5*log(Rx*Rx+Rx) )),60)  
                
            hours = wall/3600 
            wall = wall - hours * 3600
            minutes = wall/60
            wall = wall - minutes * 60
            seconds = wall
            
            if seconds == 0 or seconds == 1: 
                    seconds = 2

            hrs = str(hours)
            if len(hrs) == 1:
                hrs = '0'+hrs

            mns = str(minutes)
            if len(mns) == 1:
                mns = '0'+mns
                    
            secs = str(seconds)
            if len(secs) == 1:
                secs = '0'+secs
                    
            walltime = hrs+':'+mns+':'+secs
                
            nodes = '1'
            cores = '8'
            fl.write( '#!/bin/sh \n\n' ) 
            fl.write( '#PBS -l walltime='+walltime+'\n')
            fl.write( '#PBS -l nodes=1:ppn='+cores+'\n')
            fl.write( '#PBS -l mem='+mem+'\n')
            fl.write( '#PBS -j oe\n')
            fl.write( '#PBS -N '+jobname+'\n')
            fl.write( '#PBS -M parzuchowski@frib.msu.edu\n')
            fl.write( '#PBS -m a \n\n') 


            fl.write( 'cd $HOME/qd_imsrg/src/ \n\n')
            fl.write( 'export OMP_NUM_THREADS='+cores+'\n\n' ) 
       
            fl.write('./calc_TBME '+carg) 
            fl.write('qstat -f ${PBS_JOBID}')
            fl.close()
	
            os.system("chmod 0755 run_TBME.bat" ) 

else:
    print 'FAIL: does not compute' 
    raise SystemExit


fq.close()


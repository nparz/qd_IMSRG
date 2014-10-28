n=raw_input('number of particles: ')
hw=raw_input('hw: ')
R = raw_input('shells: ')
hw=hw.replace('.','')
hw=hw.replace('0','')

specify=n+hw+R

fmap = open('map_spec.dat','r')
fspec = open('spectrum.dat','r')
fres = open('results_'+specify+'.dat','r')
ftab = open('table_'+specify+'.dat','w') 

cor = 1.0 #for relative error
dmap = {} # dmap[Ml][Ms] = location in spectrum file 

Msmax = 0
Mlmax = 0 
for line in fmap:
    a = line.strip().split()
    a1 = int(a[0])
    a2 = int(a[1])
    a3 = int(a[2])
    if a1 > Mlmax:
        Mlmax = a1
    if a2 > Msmax:
        Msmax = a2 
        
    if a1 in dmap.keys():
        dmap[a1][a2]=a3
    else:
        dmap[a1]={}
        dmap[a1][a2]=a3

for line in fspec:
    g=line.strip().split()[1:]

resdict={} #resdict[Ml][Ms] = list of veronika levels
for line in fres:
    if ',' in line:
        a = line.strip().split(',')
        Ml = int(a[0])
        Ms = int(a[1])
        
        if Ml in resdict.keys():
            resdict[Ml][Ms]=[]
        else:
            resdict[Ml]={}
            resdict[Ml][Ms] = [] 
        continue
    if line.strip() == '':
        continue
    resdict[Ml][Ms].append(float(line.strip()))
    



for Ml in range(0,Mlmax+1) :
    for Ms in range(0,Msmax+1,2): 

        outstr = str(Ml) + ',' + str(Ms) + '\n\n'
        ftab.write(outstr) 
        
        q = resdict[Ml][Ms]
    
        if Ml==0 and Ms==0:
            e = float(g[0]) 
            b = q[0]
            per = abs(e - b)/cor

            outstr = str(e)+ '  ' + str(b) + '  '+str(per)+'\n' 
            ftab.write(outstr)
            q.remove(b)
            
        start = dmap[Ml][Ms]-2
        if Ml == Mlmax:
            if Ms+2 > Msmax:
                end = len(g) 
            else:
                end = dmap[-Mlmax][Ms+2]-2
        else:
            end = dmap[Ml+1][Ms]-2 
        
      
        

        for x in g[start:end]:
            e = float(x)
            diff = 10.
            for j in q:
                if abs(e - j) < diff:
                    diff = abs(e-j) 
                    b = j
                else:
                    break
            
            per = abs(e - b)/cor
            outstr = x+ '  ' + str(b) + '  '+str(per)+'\n' 
            ftab.write(outstr)
        
            q.remove(b)
            
        ftab.write('\n')
    
            
                
fmap.close()
ftab.close()
fres.close()
fspec.close()           
         
            


    

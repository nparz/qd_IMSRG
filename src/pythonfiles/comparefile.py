
R = 4# number of shells
A = R*(R+1) 

qn = []

# Nathan's shitty indexing method
for ms in range(-1,2,2):
    for m in range(R):    
        for n in range((R-1-m)/2+1): 
           x = []
           x.append(n)
           x.append(-1*m) 
           x.append(ms) 

           qn.append(x)
           
        if m == 0:
            continue
           
        for n in range((R-1-m)/2+1): 
           x = []
           x.append(n)
           x.append(m) 
           x.append(ms)

           qn.append(x)
       
# Veronika's alternate indexing method 
qv = []   
for E in range(1,R+1): 
    for m in range(-1*R+1,R):
        for n in range(R):
            for ms in range(-1,2,2):
            
                if E == 2*n + abs(m) + 1: 
                    x = []
                    x.append(n)
                    x.append(m) 
                    x.append(ms)

                    qv.append(x) 
                
# comparing them
mapto = []

for i in range(A):
    
    for j in range(A):
        
        if qn[j][2] == qv[i][2]:
            if qn[j][1] == qv[i][1]:
                if qn[j][0] == qv[i][0]:
                    
                    mapto.append(j+1) 
                    break

# mapto[x] = y    x is the veronika index , y is the nathan index

fei_elems = {}
feifile = open('feifile.dat','r')
	
for line in feifile:
    a = line.strip().split()
    tup = (mapto[int(a[0])],mapto[int(a[1])],mapto[int(a[2])],mapto[int(a[3])])
    fei_elems[tup] = float(a[4])
	
feifile.close()

nay_elems={}
nayfile = open('shit.dat','r')
for line in nayfile:
    a = line.strip().split()
    tup = (int(a[0]),int(a[1]),int(a[2]),int(a[3]))
    nay_elems[tup] = float(a[4])
    
nayfile.close()


for key in nay_elems.keys():
    try:
        err = abs(nay_elems[key] - fei_elems[key])
        if err > 1e-8:
            print key,nay_elems[key],fei_elems[key]
    except KeyError:
        print key, 'no have',nay_elems[key]
        g = raw_input('')

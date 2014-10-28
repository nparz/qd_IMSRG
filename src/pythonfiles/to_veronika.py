import sys

info = sys.argv 
R = int(info[1]) # number of shells
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


fl=open('elements.dat','r')
#read matrix elements from file
mat_elems = {} 
for line in fl:
    a = line.strip().split()
    
    p = int(a[0])
    q = int(a[1])
    r = int(a[2])
    s = int(a[3])
    
    g = (p,q,r,s) 
    H = a[4] 
    
    mat_elems[g] = H
    

fl.close()

#fx = open('../../veronika_CI/MasterCode-ca43a04/hamfile.dat','w')
fx = open('hamfile.dat','w')
for i in range(A):
    for j in range(i+1,A):
        for k in range(A):
            for l in range(k+1,A):
                
               p = mapto[i]
               q = mapto[j]
               r = mapto[k]
               s = mapto[l]
             
               fx.write( mat_elems[(p,q,r,s)]+'\n' ) 
               
fx.close()


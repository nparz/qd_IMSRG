import operator
fl = open('errors_eigenvecs.dat','r')

a = fl.readlines()
fl.close()

info = []
g = a[0].strip().split()
n = len(g)-2
rx = range(n+1)
ex=[]
for i in range(n+1):
    ex.append([])
errs=[]
for line in a:
    g = line.strip().split()
    for i in range(n+1):
        ex[i].append(float(g[i])) 
    errs.append(float(g[-1]))
    
x= max(errs)

# normalize errors to 1
for i in range(len(errs)): 
    errs[i] = errs[i]/x

fx = open( 'error_v_1p1h.dat', 'w')

for i in range(len(errs)): 
    min_index, min_value = min(enumerate(ex[1]), key=operator.itemgetter(1))
    perstr = str(ex[1][min_index])
    outstr = perstr+(20-len(perstr))*' '+str(errs[min_index])+'\n'
    fx.write(outstr) 
    ex[1][min_index] = 2.0
    
fx.close()
    


'''
for i in range(len(errs)):
    if (errs[i] > 0.3):
        for j in range(n+1):
            rx[j] = ex[j][i]
            
        index, value = max(enumerate(rx), key=operator.itemgetter(1))
         
        excit_lab = str(value)[:6] +' '+str(index)+'p'+str(index)+'h, '
        
        rx[index] = 0.0 
        
        index, value = max(enumerate(rx), key=operator.itemgetter(1))
         
        excit_lab += str(value)[:6] +' '+str(index)+'p'+str(index)+'h, '
        
        rx[index] = 0.0
        
        index, value = max(enumerate(rx), key=operator.itemgetter(1))
        
        excit_lab += str(value)[:6] +' '+str(index)+'p'+str(index)+'h'
        
        outstr = 'error: '+str(errs[i])[:3] +' '
        outstr += '    composition: '+excit_lab
        print outstr
        
            
'''



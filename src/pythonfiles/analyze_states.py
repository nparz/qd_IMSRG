fl = open('level_percentages.dat','r')
fj = open('CI_spectrum.dat','r') 

a = fl.readlines()

g = a[0].strip().split() 

n = len(g)-1

perc_info = [] 

for line in a: 

    g = line.strip().split()
    q = [] 
    for element in g:
        q.append(float(element)) 
    perc_info.append(q) 
   
fl.close()

CI_info=[]
for line in fj:
    g = line.strip().split()
    q = []
    for element in g:
        q.append(float(element))
    CI_info.append(q)
  
fj.close()

fx = open('errors_eigenvecs.dat','w') 
R = len(perc_info) 

for i in range(R): 
    error = str(abs(CI_info[0][i+1] - CI_info[-1][i+1])/CI_info[0][1])
    m = len(error)
    if 'e' in error:
        error = error[:6]+error[-4:] 
    elif (m > 10):
        error = error[:10]
    elif (m < 10): 
        error = error + (10-m)*'0'
    outstr = ''
    for el in perc_info[i]:
        kx = str(el) 
        m = len(kx) 
        if 'e' in kx:
            kx = kx[:6]+kx[-4:] 
        elif (m > 10):
            kx = kx[:10]
        elif (m < 10): 
            kx = kx + (10-m)*'0'
        outstr += kx + '  ' 
    outstr += error +'\n'

    fx.write(outstr)

fx.close()

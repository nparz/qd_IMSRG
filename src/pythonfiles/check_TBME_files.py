fl=open('cases.txt','r')
flag='off'

for line in fl:
    a=line.split()
    for i in range(3):
        a[i]=a[i].strip()
    
    b=a[0][0:4]
    b=b.replace('d','0')
    if b[0]=='.':
        b='0'+b[0:3]
    fn='../TBMEfiles/TBME_'+b+'_'+a[1]+'_'+a[2]+'.dat'

    try: 
        with open(fn): pass
    except IOError:
        if flag == 'off':
            fx=open('missingTBMEfiles.txt','w') 
            flag='on'
        print 'file: '+fn[13:]+' absent.\n'
        print 'submitting for TBME calculation'
        fx.write(line)        
    
fl.close()
if flag == 'on':
    fx.close()
    

    

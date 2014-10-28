fl = open('flowing_CI.dat','r')
fx = open('cut_flow_CI.dat','w')

for line in fl:
    a = line.strip().split()
    outline = ''
    for i in range(11):
        outline = outline + a[i] + ' '
    outline += '\n'
    fx.write(outline)

fx.close()
fl.close()

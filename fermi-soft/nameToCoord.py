def run(File):
    f = open(File,'r')
    lines = f.readlines()
    result=[]
    for x in lines:
        result.append(x.split(',')[0])
    f.close()

    i=0
    while i < len(result):
        W = result[i]
        Rh = float(W[0:2])
        Rm = float(W[2:4])
        Rs = float(W[4:8])
        sign = W[8]
        Dd = float(W[9:11])
        Damin = float(W[11:13])
        Dasec = float(W[13:15])
        RA = (Rh * 15) + (Rm * 0.25) + (Rs * (15.0/3600))
        dec = Dd + (Damin * 1.0/60) + (Dasec * 1.0/3600)
        if sign == '+':
            DEC = dec
        if sign == '-':
            DEC = -1 * dec
        print str(RA) + "\t" + str(DEC)
        i = i + 1

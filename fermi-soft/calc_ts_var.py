import numpy as np
import matplotlib.pyplot as plt

def run(Fname):

    if Fname[0] == "_":
        f = open("/storage/home/mwt5345/work/Analysis/2WHSP/sources/all.txt")
        f = f.readlines()
        i = 0
        while i < len(f):
            if "_" + Fname.split("_")[1] == f[i].split("\t")[1].split("\n")[0]:
                Name = f[i].split("\t")[0]
            i += 1
    else:
        Name = Fname.split("_")[0]
   
    F = open("../results.txt")
    F = F.readlines()

    flu = 0

    i = 0
    while i < len(F):
        holdi = F[i].split("\t")
        if Name == holdi[0]:
            flu = float(holdi[2])
        i += 1

    T = np.genfromtxt(Fname,delimiter=',',usecols=7)#
    UL = np.genfromtxt(Fname,delimiter=',',usecols=4)#

    dE = np.genfromtxt(Fname,delimiter=',',usecols=1)#
    E = np.genfromtxt(Fname,delimiter=',',usecols=0)#
    TestS = np.genfromtxt(Fname,delimiter=',',usecols=5)#
    LogL = np.genfromtxt(Fname,delimiter=',',usecols=9)#
    LogL0 = np.genfromtxt(Fname,delimiter=',',usecols=10)#

    Tul=[]
    Tf=[]
    Ulul=[]
    dEf=[]
    Ef=[]
    dEul=[]
    Eul=[]
    LLul=[]
    LL0ul=[]
    LLf=[]
    LL0f=[]

    i = 0
    while i < len(TestS):
        if TestS[i] < 7:
            #Tul.append(MJD[i])
            Ulul.append(UL[i])
            dEul.append(dE[i])
            Eul.append(E[i])
            LLul.append(LogL[i])
            LL0ul.append(LogL0[i])
        else:
            #Tf.append(MJD[i])
            Ef.append(E[i])
            dEf.append(dE[i])
            LLf.append(LogL[i])
            LL0f.append(LogL0[i])
        i = i + 1

    #Tul = np.array(Tul)
    #Tf = np.array(Tf)
    Ulul = np.array(Ulul)
    dEf = np.array(dEf)
    Ef = np.array(Ef)
    LLf = np.array(LLf)
    dEul = np.array(dEul)
    Eul = np.array(Eul)
    LLul = np.array(LLul)
    LL0ul = np.array(LL0ul)
    LL0f = np.array(LL0f)

    # Make sure to updata TS var to incude more than just ULs
    #TSvar = abs(2 * (TotTS - avTS))
    tsn = LLf - LL0f
    stuff = ((dEf)**2)/((dEf)**2 + (0.02**2)*(flu**2))
    TSvar1 = 2 * np.sum(np.power((stuff),1)*tsn)


    tsm = LLul - LL0ul
    stuff = (0.5*(Ulul - Eul))**2/((0.5*(Ulul - Eul))**2 + (0.02**2)*(flu**2))
    TSvar2 = 2 * np.sum(np.power((stuff),1)*tsm)

    TSvar = abs(TSvar1 + TSvar2)
    return TSvar

f = open("files.txt")
f = f.readlines()
print len(f)
i = 0
hold = np.zeros(len(f))
while i < len(f):
    v = f[i].split("\n")[0]
    t = run(v)
    hold[i] = t
    print str(v.split("_")[0]) + "\t" + str(t)
    i += 1


plt.hist(hold,bins=100,normed=False)

plt.title(r'TS$_{var}$ Distribution')
plt.xlabel(r'TS$_{var}$')
plt.ylabel('Number of Sources')
plt.axvline(72.44,c='r',lw=2)
plt.show()

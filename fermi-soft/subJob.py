from xml.etree import ElementTree as ET
import subprocess
import numpy as np

def make(name,NN):
    mod = ET.parse(str(name) + '_output_model.xml')
    model = mod.getroot()
    source = model.getchildren()

    #Find the source with the correct name
    i = 0
    while i < len(source):
        if model[i].attrib['name'] == str(name):
                params = source[i].getchildren()
        else:
            pass
        i += 1

    coords = params[1]
    RA = coords[0].attrib['value']
    DEC = coords[1].attrib['value']

    xx = "#!/bin/bash\n#PBS -l walltime=24:00:00\n#PBS -l nodes=1:ppn=10\n\ncd ~/scratch/Analysis/2WHSP/JOBS/\n\npython -c \'import multiMakeUBLikeLC as m;m.run(\"" + str(name) + "\","+str(RA)+","+str(DEC)+",1e3,3e5,\"../spacecraft/spacecraft/spacecraft.fits\",8,0.1,239587201,496426332.0,\"@../extended/extended/photons.txt\",1,100,128,3,1e6,0.01,50,False,10)\'"

    f = open("Z" + str(NN) + "_lc.sh","w")
    f.write(xx)

    #subprocess.call(['qsub',str(N) + '_lc.sh'])

xx = np.genfromtxt('/storage/home/mwt5345/work/Analysis/2WHSP/sources/all.txt',delimiter='\t',dtype=str,usecols=1)

i = 0

while i < len(xx):
    make(xx[i],i)
    i = i + 1

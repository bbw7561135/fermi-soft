"""

    runMany.py - This script does a complete Fermi unbinned likelihood
                    analysis. Currently the script is configured to
                    handle the 2WHSP catalog data.

    Author - Michael W. Toomey
    Institution - The Pennsylvania State University

"""

import subprocess
import makeUBFermiFiles2 as mF
import multiUBLike2 as mL
from make3FGLxml import *
import multiprocessing
import numpy as np

def runMany(num,Name):
    X = np.genfromtxt('/storage/home/mwt5345/work/Analysis/2WHSP/sources/' + str(Name),names=True,delimiter='\t',dtype=None)#,skiprows=2)
    i = 0
    pool = multiprocessing.Pool(num)
    while i < len(X):
        Y = X[i]
        if False:
            pass
        else:
            if False:
                pass
            else:
                name = Y[0]
                FGL = Y[1]
                pool.apply_async(run,(name,FGL,))
        i = i + 1
    pool.close()
    pool.join()
    print 'Data has been compiled. Please see results.txt for results.'

"""Program that runs entire Unbinned Likelihood analysis without interuption.
Optimized for detection of faint sources at high energy, > 1 GeV."""
def run(Name,FGL,minEnergy=1e3,maxEnergy=3e5,SCFile="../spacecraft/spacecraft/spacecraft.fits",radius=8,binsz=0.1,TSTART=239587201,TSTOP=496426332.0,Infile="@../extended/extended/photons.txt",bins=1,zmax=100,evclass=128,evtype=3,ulTS=25,NMtol=0.01):
    WHSPname = Name
    coords = nametocoord(Name)
    if '3FGL' in FGL:
        Name = FGL
    else:
        pass
    #makeUBFermiFiles
    mF.run(Name,coords[0],coords[1],minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,TSTOP,Infile,bins,zmax,evclass,evtype)
    #XML model creation
    makeRawXML(Name)
    #adding source of interest to xml model
    if '3FGL' in FGL:
        pass
    else:
        addSourceXML(Name,coords[0],coords[1])
    #multiUBLike - returns an array that looks like this
    #               Return = [log-likelihood,flux,flux error,Baeysian upper limit,test statisitc, and others]
    fit = mL.run(Name,coords[0],coords[1],minEnergy,maxEnergy,SCFile,ulTS,NMtol,evclass,str(Name)+'_in_model.xml')
    printToFile(WHSPname,fit)

"""
    For sources that are not included in the generic 3FGL catologue, this script
    adds a new model to the file for a source with a simple power law.
"""
def addSourceXML(Name,Ra,Dec):
    #For sources that are not included in the generic 3FGL catologue, this script
    #adds a new model to the file for a source with a simple power law.
    f = open(str(Name)+"_in_model.xml","r")
    y = f.readlines()
    f.close()
    string = "<source ROI_Center_Distance=\"0.00\" name=\"" + str(Name) + "\" type=\"PointSource\">\n\
\t<spectrum apply_edisp=\"false\" type=\"PowerLaw\">\n\
\t\t<parameter free=\"1\" max=\"1e7\" min=\"1e-7\" name=\"Prefactor\" scale=\"1e-12\" value=\"1.0\"/>\n\
\t\t<parameter free=\"1\" max=\"10.0\" min=\"0.0\" name=\"Index\" scale=\"-1.0\" value=\"2.0\"/>\n\
\t\t<parameter free=\"0\" max=\"5e5\" min=\"30\" name=\"Scale\" scale=\"1.0\" value=\"3000\"/>\n\
\t</spectrum>\n\
\t<spatialModel type=\"SkyDirFunction\">\n\
\t\t<parameter free=\"0\" max=\"360.0\" min=\"-360.0\" name=\"RA\" scale=\"1.0\" value=\"" + str(Ra) + "\"/>\n\
\t\t<parameter free=\"0\" max=\"90\" min=\"-90\" name=\"DEC\" scale=\"1.0\" value=\"" + str(Dec) + "\"/>\n\
\t</spatialModel>\n\
</source>"
    y.insert(4,string)
    f = open(str(Name)+"_in_model.xml","w")
    y="".join(y)
    f.write(y)
    f.close()

"""
    Makes a generic XML model for the region of interest using make3FGLxml
"""
def makeRawXML(Name):
    mymodel = srcList('/storage/home/mwt5345/work/Fermi/Z-RandomShit/AnalysisFiles/gll_psc_v16.fit', Name + '_gtmktime.fits', str(Name) + '_in_model.xml')
    mymodel.makeModel(normsOnly=False,extDir='',radLim=3,maxRad=5,sigFree=4,varFree=True,oldNames=True,ExtraRad=5)

"""
    Prints the result of the likelihood to a file
"""
def printToFile(Name,arr):
    f = open('results.txt','a')
    f.write(Name + '\t')
    for i in arr:
        f.write("%s\t" % i)
    f.write('\n')
    f.close()

"""
    Converts the name of the source to Ra and Dec in degrees
"""
def nametocoord(Name):
    W = Name
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
    return [RA,DEC]

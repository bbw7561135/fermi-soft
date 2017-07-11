import sys
sys.path.insert(0,'/Projects/GLAST/share/bin')
sys.path.insert(0,'/home/mtoomey/FilesAndScripts')
from FermiObject import FermiObject
from BinnedAnalysis import *
import numpy as np
import subprocess
from math import sqrt

'''

    Name    # The name of your source <String>
    Ra      # The right ascension of your source <int/float>
    Dec     # The declination of your source <int/float>

'''

def runMakeResMap(Name,Ra,Dec,minEnergy,maxEnergy,radius,binsz):
    print "This is makeResMap.\nStarting construction of residuals map.\n"
    f = FermiObject()
    print "Creating model map of the region surrounding " + str(Name) + " with gtmodel."
    f._setSrcmaps(Name + '_srcmaps.fits')
    f._setExpcube( '../' + Name + '_ltcube.fits')
    f._setBexpmap( '../' + Name + '_expcube.fits')
    f._setSrcmdl(Name + '_output_model.xml')
    f._setOutfile( Name + '_model_map.fits')
    f._setIrfs('CALDB')
    f.amonModelmaps()
    print "Finished model map. Now constructing counts map for comparison."
    f._setAlgorithm('CMAP')
    f._setEvfile( '../' + Name + '_gtmktime.fits')
    f._setOutfile( Name + '_cmap_small.fits')
    f._setScfile('NONE')
    pix = int((sqrt(2)*radius)/float(binsz))
    f._setNxpix(pix)
    f._setNypix(pix)
    f._setBinsz(binsz)
    f._setCoordsys('CEL')
    f._setAxisrot(0)
    f._setProj('AIT')
    f._setRa(Ra)
    f._setDec(Dec)
    f._setEmin(minEnergy)
    f._setEmax(maxEnergy)
    f.amonBincmap()
    print "Now making residuals map with ftool farith."
    yes = subprocess.call(['farith', Name + '_cmap_small.fits',Name + '_model_map.fits',Name + '_residual.fits','SUB'])
    if yes == 0:
        print "Residuals map has been created, see " + str(Name) + "_residuals.fits"
    else:
        print "Subprocessing failed. Unable to create residuals map with farith."

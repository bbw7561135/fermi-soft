"""

	makeUBFermiFiles.py - client code for FermiObject class that will do
	all the steps required for an unbinned analysis upto the creation of the .xml 
    model. This is done becasue the model should be looked over by human eyes 
    before continuing forward. In addition, sources you wish to keep fixed or
    free should be carefully concidered before proceeding forward with a 
    likelihood analysis.

	Author: Michael Toomey <mwt5345@psu.edu>
	Institution: The Pennsylvania State University
    
    Parameters:
    -------------------------------------------------------------------------------
    Name           # The name of your source <String>
    RA             # The right ascension of your source in degrees <float,int>
    DEC            # The declination of your source in degrees <float,int>
    minEnergy      # The minimum energy of your analysis in MeV <float,int>
    maxEnergy      # The maximum energy of your analysis in MeV <float,int>
    SCFile         # The spacecraft file with full path for the analysis <String>
    radius         # The radius of your region of interest in degrees <float,int>
    binsz          # The size of your pixels for your various maps in degrees
                        (recommend keeping 0.1 degrees/pixel) <float,int> 
    TSTART         # The starting time of your analysis in MET <float,int>
    TSTOP          # The ending time of your analysis in MET <float,int>
    Infile         # The raw photon file with full path for the analysis <String>
    bins           # The number of jobs you want to split the analysis into while 
                     making the livetime <int>
    zmax           # The maximum apparent zenith angle in degree <float,int>
    evclass        # The event class for the analysis <int>
    evtype         # The event type for the analysis <int>

"""

import sys
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi/gtapps_mp-1.1')
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi')
from gtltcube_mp import *
from FermiObject import FermiObject
from math import log10,sqrt
import gt_apps as my_apps

def run(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,TSTOP,Infile,bins,zmax,evclass,evtype):
    
    print "This is makeFermiFiles."   

    f = FermiObject()

    if evtype == 512:
        irf = "P8R2_ULTRACLEAN_V6"
    elif evtype == 128:
        irf = "P8R2_SOURCE_V6"
    elif evtype == 256:
        irf = "P8R2_CLEAN_V6"
    elif evtype == 1024:
        irf = "P8R2_ULTRACLEANVETO_V6"

    """

	    Following steps execute Fermi Tool gtselect

    """

    print('\nWorking on file.')
    print('Cutting file to fit desired parameters . . .\n')
    f._setEvclass(evclass)
    f._setEvtype(evtype)
    f._setRa(RA)
    f._setDec(DEC)
    f._setRad(radius)
    f._setEmin(minEnergy)
    f._setEmax(maxEnergy)
    f._setZmax(zmax)
    f._setTmin(TSTART)	
    f._setTmax(TSTOP)
    f._setInfile(Infile)
    f._setOutfile( Name + '_gtselect.fits')
    f.amonSelect()
    print('File cuts have been made. Now making cuts for GTI using spacecraft file.')

    """

	    Following steps execute Fermi Tool gtmktime

    """

    f._setScfile(SCFile)
    f._setRoicut('no')
    f._setEvfile( Name + '_gtselect.fits')
    f._setOutfile( Name + '_gtmktime.fits')
    ###############################################
    #         Filter expression                   #
    Filter = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
    ###############################################
    f._setFilter(Filter)
    print('Working on file ' + str(f.getOutfile()) + '. . .')
    f.amonTime()
    print('File cuts have been made. Now begining construction of the counts map from event data.')

    """

	    Following steps execute Fermi Tool gtbin
	    to create the counts map

    """

    f._setAlgorithm('CMAP')
    f._setEvfile( Name + '_gtmktime.fits')
    f._setOutfile( Name + '_cmap.fits')
    f._setScfile('NONE')
    num_pix = int((2*radius)/float(binsz))
    print "Counts map is " + str(num_pix) + " by " + str(num_pix) + " pixels."
    f._setNxpix(num_pix)
    f._setNypix(num_pix)
    f._setBinsz(binsz)
    f._setCoordsys('CEL')
    f._setAxisrot(0)
    f._setProj('AIT')
    f.amonBincmap()
    print('Counts map is complete. Now begining construction of the livetime cube.')

    #You will want to run gtltcube_mp.py intead of the standard
    #Fermi Tool gtltcube as this script enables multiprocessing
    #capabilities that greatly increases the speed.
    print "Now working on ltcube file using gtltcube_mp.py.\n"
    gtltcube_mp(bins, SCFile, Name + '_gtmktime.fits', Name + '_ltcube.fits', False, zmax)
    print "\nltcube complete.\nMoving to compute exposure map with gtexpmap.\n"

    print "\nltcube complete.\nMoving to compute exposure map with gtexpmap.\n"
    my_apps.expMap['evfile'] = Name + '_gtmktime.fits'
    my_apps.expMap['scfile'] = SCFile
    my_apps.expMap['expcube'] = Name + '_ltcube.fits'
    my_apps.expMap['outfile'] = Name + '_expMap.fits'
    my_apps.expMap['irfs'] = 'CALDB'
    my_apps.expMap['srcrad'] = radius + 10
    my_apps.expMap['nlong'] = 4*(radius + 10)
    my_apps.expMap['nlat'] = 4*(radius + 10)
    ebin = int(10*log10(maxEnergy/minEnergy))
    print "There are " + str(ebin) + " energy bans."
    my_apps.expMap['nenergies'] = ebin
    my_apps.expMap.run()
    print "Finnished making exposure map.\n"

def cli():
    import argparse
    helpString="\nParameters - \
    Name: The name of your source <String>,\
    RA: The right ascension of your source in degrees <float,int>,\
    DEC: The declination of your source in degrees <float,int>,\
    TSTART: The starting time of your analysis in MET <float,int>,\
    TSTOP: The ending time of your analysis in MET <float,int>\n \
    infile: The raw photon file with full path for the analysis <String>,\
    bins: The number of jobs you want to split the analysis into while\
    making the livetime <int>,\
    radius: The radius of your region of interest in degrees <float,int>,\
    minEnergy: The minimum energy of your analysis in MeV <float,int>,\
    maxEnergy: The maximum energy of your analysis in MeV <float,int>,\
    binsz: The size of your pixels for your various maps in degrees\
    (recommend keeping 0.1 degrees/pixel) <float,int>,\
    zmax: The maximum apparent zenith angle in degrees <float,int>,\
    SCFile: The spacecraft file with full path for the analysis <String>"
    parser=argparse.ArgumentParser(description=helpString)
    args=parser.parse_args()

if __name__=='__main__': cli()

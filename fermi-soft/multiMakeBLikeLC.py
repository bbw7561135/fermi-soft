"""

        makeBLikeLC.py - creates a light curve for a source of interest
        using a binned likelihood analysis. Be sure to have a spectral model for
        your source of interest in your output model from your likelihood analysis
        for the entire dataset; you will be using this XML file throughout.

        Author: Michael Toomey <mwt5345@psu.edu>
        Institution: The Pennsylvania State University

        * Relies on the implementation of the FermiObject.py to gain access to the
            use of the various Fermi Tools required for this analysis.

"""

from FermiObject import FermiObject
import sys
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi/gtapps_mp-1.1') #For use at NRL replace with path to directory with gtltcube_mp
from gtltcube_mp import *
from BinnedAnalysis import *
import IntegralUpperLimit as IUL
import numpy as np
import subprocess
from math import sqrt,log10
import gt_apps as my_apps
import pyfits
import multiprocessing
from ConfigParser import ConfigParser

# Conducts the complete likelihood analysis to get the flux, flux error, and
# TS value for each bin of the light curve as well as upper limits when needed.
def runFermiTools(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,TSTOP,Evfile,bins,zmax,evclass,evtype,TSul,NMtol,lc_bin_num,runMRM):

    print "Working on bin " + str(lc_bin_num) + " for the light curve."

    f = FermiObject()

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
    f._setInfile(Evfile)
    f._setOutfile( Name + '_gtselect' + str(lc_bin_num) + '_lc.fits')
    f.amonSelect()
    print('File cuts have been made. Now making cuts for GTI using spacecraft file.')

    """

        Following steps execute Fermi Tool gtmktime

    """

    f._setScfile(SCFile)
    f._setRoicut('no')
    f._setEvfile( Name + '_gtselect' + str(lc_bin_num) + '_lc.fits')
    f._setOutfile( Name + '_gtmktime' + str(lc_bin_num) + '_lc.fits')
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
    f._setEvfile( Name + '_gtmktime' + str(lc_bin_num) + '_lc.fits')
    f._setOutfile( Name + '_cmap' + str(lc_bin_num) + '_lc.fits')
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
    print('Counts map is complete.Now begining construction of the counts cube.')



    """

        Following steps execute Fermi Tool gtbin
        to create counts cube (3D counts map).

    """

    f._setAlgorithm('CCUBE')
    f._setOutfile( Name + '_ccube' + str(lc_bin_num) + '_lc.fits')
    pix = int((sqrt(2)*radius)/float(binsz))
    print "Counts cube is " + str(pix) + " by " + str(pix) + " pixels."
    f._setNxpix(pix)
    f._setNypix(pix)
    ebin = int(10*log10(maxEnergy/minEnergy))
    print "There are " + str(ebin) + " logarithmically uniform energy bins."
    f._setEnumbins(ebin)
    f.amonBinccube()
    print('Counts cube is complete.')

    print('Using XML model from whole dataset.\n Moving on to gtltcube.')

    print "Now working on ltcube file using gtltcube\n"

    my_apps.expCube['evfile'] = Name + '_gtmktime' + str(lc_bin_num) + '_lc.fits'
    my_apps.expCube['scfile'] =  SCFile
    my_apps.expCube['outfile'] = Name + '_ltcube' + str(lc_bin_num) + '_lc.fits'
    my_apps.expCube['dcostheta'] = 0.025
    my_apps.expCube['binsz'] = 1
    my_apps.expCube['phibins'] = 0
    my_apps.expCube['zmax'] = zmax
    my_apps.expCube['chatter'] = 0
    my_apps.expCube.run()

    print "\nltcube complete.\nMoving to compute exposure map for whole sky with gtexpcube2.\n"
    
    f._setInfile(Name + '_ltcube' + str(lc_bin_num) + '_lc.fits')
    cubePix = int((2*radius + 20)/binsz)
    f._setNxpix(cubePix)
    f._setNypix(cubePix)
    f._setBinsz(binsz)
    f._setCoordsys('CEL')
    f._setRa(RA)
    f._setDec(DEC)
    f._setAxisrot(0)
    f._setProj('AIT')
    f._setEmin(minEnergy)
    f._setEmax(maxEnergy)
    f._setEnumbins(ebin)
    f._setOutfile( Name + '_expcube' + str(lc_bin_num) + '_lc.fits')
    f._setIrfs('P8R2_SOURCE_V6')
    f.amonExpcube2()

    print "Finnished making exposure map.\n"

    print "Now begining computation of model counts map with gtsrcmaps.\n"

    f._setScfile(SCFile)
    f._setExpcube( Name + '_ltcube' + str(lc_bin_num) + '_lc.fits')
    f._setCmap( Name + '_ccube' + str(lc_bin_num) + '_lc.fits')
    f._setBexpmap( Name + '_expcube' + str(lc_bin_num) + '_lc.fits')
    f._setOutfile( Name + '_srcmaps' + str(lc_bin_num) + '_lc.fits')
    f._setIrfs('CALDB')
    f._setSrcmdl(Name + '_output_model.xml')
    f.amonSrcmaps()

    print "Finished srcmaps. Now moving to conduct a BINNED likelihood analysis."

    obs = BinnedObs(binnedExpMap= Name + '_expcube' + str(lc_bin_num) + '_lc.fits',expCube= Name + '_ltcube' + str(lc_bin_num) + '_lc.fits',srcMaps= Name + '_srcmaps' + str(lc_bin_num) + '_lc.fits', irfs='P8R2_SOURCE_V6')
    analysis = BinnedAnalysis(obs,srcModel= Name + '_output_model.xml', optimizer='NEWMINUIT')
    likeObj = pyLike.NewMinuit(analysis.logLike)
    analysis.tol = NMtol
    analysis.fit(verbosity=0,covar=True,optObject=likeObj)
    fit = likeObj.getRetCode()
    print "Likelihood has converged whith Code " + str(likeObj.getRetCode())
    Flux = analysis.flux(Name, emin=minEnergy, emax=maxEnergy)
    Ferr = analysis.fluxError(Name, emin=minEnergy, emax=maxEnergy)
    MeVtoErg = 1.602e-6
    ef = analysis.energyFlux(Name,minEnergy,maxEnergy)*MeVtoErg
    ef_err =  analysis.energyFluxError(Name,minEnergy,maxEnergy)*MeVtoErg
    UL = False
    TSUM = TSTART + TSTOP
    TMID = TSUM/2
    limit = Flux
    if analysis.Ts(Name) < TSul:
        UL = True
        limit,results = IUL.calc_int(analysis, Name,cl=0.90, emin=minEnergy, emax=maxEnergy)

    #Run gtselect to make smaller data fits file to compute the exposure, set to 3 degrees around source of interest
    f._setRad(3)
    f._setInfile(Evfile)
    f._setOutfile( Name + '_gtselect' + str(lc_bin_num) + '_exposure.fits')
    print "Creating file " + Name + "_gtselect_exposure.fits"
    f.amonSelect()

    #Run gtmaketime on this small region
    f._setEvfile( Name + '_gtselect' + str(lc_bin_num) + '_exposure.fits')
    f._setOutfile( Name + '_gtmktime' + str(lc_bin_num) + '_exposure.fits')
    print('Working on file ' + str(f.getOutfile()))
    f.amonTime()

    my_apps.evtbin['algorithm'] = 'LC'
    my_apps.evtbin['evfile'] = f.getOutfile()
    my_apps.evtbin['outfile'] = Name + '_LC' + str(lc_bin_num) + '_exposure.fits'
    my_apps.evtbin['scfile'] = f.getScfile()
    my_apps.evtbin['tbinalg'] = 'LIN'
    my_apps.evtbin['tstart'] = f.getTmin()
    my_apps.evtbin['tstop'] = f.getTmax()
    my_apps.evtbin['dtime'] = TSTOP - TSTART
    my_apps.evtbin.run()

    yes = subprocess.call(['gtexposure', Name + '_LC' + str(lc_bin_num) + '_exposure.fits',f.getScfile(),'P8R2_SOURCE_V6',Name + '_output_model.xml',Name])
    if yes == 0:
        print "Exposure map has been created"
    else:
        print "Subprocessing failed. Unable to create exposure map with gtexposure."
    
    print "Time bin complete."

    hdulist = pyfits.open(Name + '_LC' + str(lc_bin_num) + '_exposure.fits')
    tbdata = hdulist[1].data
    z = tbdata['EXPOSURE']
    exp = z[0]

    ################################################################
    #           This portion prints to the text file               #
    ################################################################

    f = open("lc_output.txt", "a")
    f.write(str(Flux) + ',' + str(Ferr) + ',' + str(ef) + ',' + str(ef_err) + ',' + str(limit) + ','+ str(analysis.Ts(Name)) + ',' + str(UL) + ',' + str(TMID) + ',' + str(exp) + '\n')
    f.close()
    print "Likelihood analysis on this band is complete."
    if runMRM == True:
        runMakeResMap(Name, RA, DEC, radius, binsz,minEnergy,maxEnergy,lc_bin_num)
    else:
        pass

    yes = subprocess.call(['rm',Name + '_gtselect' + str(lc_bin_num) + '_lc.fits',Name + '_gtmktime' + str(lc_bin_num) + '_lc.fits',Name + '_cmap' + str(lc_bin_num) + '_lc.fits',Name + '_ccube' + str(lc_bin_num) + '_lc.fits',Name + '_ltcube' + str(lc_bin_num) + '_lc.fits',Name + '_expcube' + str(lc_bin_num) + '_lc.fits',Name + '_LC' + str(lc_bin_num) + '_exposure.fits',Name + '_srcmaps' + str(lc_bin_num) + '_lc.fits',Name + '_gtselect' + str(lc_bin_num) + '_exposure.fits',Name + '_gtmktime' + str(lc_bin_num) + '_exposure.fits'])
    if yes == 0:
        print 'Files for bin have been deleted'
    else:
        print "Subprocessing failed. Unable to delete files for bin."

"""

        This portion of the code will split the data into bands and run a likelihood analysis
        from start to finish on the data.

"""
def run():

    parser = ConfigParser()
    parser.read('../../config.ini')

    #parameters from config file
    Name = str(parser.get('params','Name'))
    RA = float(parser.get('params','RA'))
    DEC = float(parser.get('params','DEC'))
    minEnergy = float(parser.get('params','minEnergy'))
    maxEnergy = float(parser.get('params','maxEnergy'))
    SCFile = str(parser.get('params','SCFile'))
    radius = float(parser.get('params','radius'))
    binsz = float(parser.get('params','binsz'))
    TSul = float(parser.get('params','ulTS'))
    BUL = bool(parser.get('params','BUL'))
    SRC = bool(parser.get('params','SRC'))
    Mtol = float(parser.get('params','Mtol'))
    NMtol = float(parser.get('params','NMtol'))
    runMRM = False #bool(parser.get('params','makeRes'))
    TSTART = float(parser.get('params','TSTART'))
    TSTOP = float(parser.get('params','TSTOP'))
    Evfile = str(parser.get('params','Infile'))
    bins = int(parser.get('params','bins'))
    zmax = float(parser.get('params','zmax'))
    evclass = int(parser.get('params','evclass'))
    evtype = int(parser.get('params','evtype'))
    bin_num = int(parser.get('params','bin_num'))
    num = int(parser.get('params','jobs'))

    TTime = TSTOP - TSTART
    bin_size = TTime/bin_num
    i = 0
    pool = multiprocessing.Pool(num)
    while i < bin_num:
        tstop = TSTART + bin_size
        pool.apply_async(runFermiTools, args=(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,tstop,Evfile,bins,zmax,evclass,evtype,TSul,NMtol, i + 1,runMRM))
        i = i + 1
        TSTART = TSTART + bin_size
    pool.close()
    pool.join()
    print 'LC data has been compiled. Please see lc_output.txt for results.'
"""
    Copy of makeResMap - needed to copy to avoid error from farith everytime residuals map was made
"""
def runMakeResMap(Name, Ra, Dec, radius, binsz,minEnergy,maxEnergy,lc_bin_num):
    print "This is makeResMap.\nStarting construction of residuals map.\n"
    f = FermiObject()
    print "Creating model map of the region surrounding " + str(Name) + " with gtmodel."
    f._setSrcmaps(Name + '_srcmaps' + str(lc_bin_num) + '_lc.fits')
    f._setExpcube(Name + '_ltcube' + str(lc_bin_num) + '_lc.fits')
    f._setBexpmap(Name + '_expcube' + str(lc_bin_num) + '_lc.fits')
    f._setSrcmdl(Name + '_output_model.xml')
    f._setOutfile( Name + '_model_map_' + str(lc_bin_num) + '.fits')
    f._setIrfs('CALDB')
    f.amonModelmaps()
    print "Finished model map. Now constructing counts map for comparison."
    f._setAlgorithm('CMAP')
    f._setEvfile(Name + '_gtmktime' + str(lc_bin_num) + '_lc.fits')
    f._setOutfile( Name + '_cmap_small_' + str(lc_bin_num) + '.fits')
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
    yes = subprocess.call(['farith', Name + '_cmap_small_' + str(lc_bin_num) + '.fits',Name + '_model_map_' + str(lc_bin_num) + '.fits',Name + '_residual_' + str(lc_bin_num) + '.fits','SUB'])
    if yes == 0:
        print 'Residuals map has been created, see ' + str(Name) + '_residuals_' + str(lc_bin_num) + '.fits'
    else:
        print "Subprocessing failed. Unable to create residuals map with farith."

"""

        makeBLikeLC.py - creates a light curve for a source of interest
        using a binned likelihood analysis. Be sure to have a spectral model for
        your source of interest in your output model from your likelihood analysis
        for the entire dataset; you will be using this XML file throughout.

        Author: Michael Toomey <mwt5345@psu.edu>
        Institution: The Pennsylvania State University

        Modified for use on NRL's HESELIN system

        * Relies on the implementation of the FermiObject.py to gain access to the
            use of the various Fermi Tools required for this analysis.

"""

from FermiObject import FermiObject
import sys
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi/gtapps_mp-1.1')
from gtltcube_mp import *
from BinnedAnalysis import *
import IntegralUpperLimit as IUL
import numpy as np
import subprocess
from math import sqrt,log10
import gt_apps as my_apps
import pyfits

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
    f._setOutfile( Name + '_gtselect_lc.fits')
    f.amonSelect()
    print('File cuts have been made. Now making cuts for GTI using spacecraft file.')

    """

        Following steps execute Fermi Tool gtmktime

    """

    f._setScfile(SCFile)
    f._setRoicut('no')
    f._setEvfile( Name + '_gtselect_lc.fits')
    f._setOutfile( Name + '_gtmktime_lc.fits')
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
    f._setEvfile( Name + '_gtmktime_lc.fits')
    f._setOutfile( Name + '_cmap_lc.fits')
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
    f._setOutfile( Name + '_ccube_lc.fits')
    pix = int((sqrt(2)*radius)/float(binsz))
    print "Counts cube is " + str(pix) + " by " + str(pix) + " pixels."
    f._setNxpix(pix)
    f._setNypix(pix)
    ebin = int(10*log10(maxEnergy/minEnergy))
    print "There are " + str(ebin) + " logarithmically uniform energy bins."
    f._setEnumbins(ebin)
    f.amonBinccube()
    print('Counts cube is complete.')

    print('Using XML model from whole dataset.\n Moving on to multiprocessing version of gtltcube.')

    #You will want to run gtltcube_mp.py intead of the standard
    #Fermi Tool gtltcube as this script enables multiprocessing
    #capabilities that greatly increases the speed.
    print "Now working on ltcube file using gtltcube_mp.py.\n"
    gtltcube_mp(bins, SCFile, Name + '_gtmktime_lc.fits', Name + '_ltcube_lc.fits', False, zmax)
    print "\nltcube complete.\nMoving to compute exposure map for whole sky with gtexpcube2.\n"
    
    f._setInfile(Name + '_ltcube_lc.fits')
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
    f._setOutfile( Name + '_expcube_lc.fits')
    f._setIrfs('P8R2_SOURCE_V6')
    f.amonExpcube2()

    print "Finnished making exposure map.\n"

    print "Now begining computation of model counts map with gtsrcmaps.\n"

    f._setScfile(SCFile)
    f._setExpcube( Name + '_ltcube_lc.fits')
    f._setCmap( Name + '_ccube_lc.fits')
    f._setBexpmap( Name + '_expcube_lc.fits')
    f._setOutfile( Name + '_srcmaps_lc.fits')
    f._setIrfs('CALDB')
    f._setSrcmdl(Name + '_output_model.xml')
    f.amonSrcmaps()

    print "Finished srcmaps. Now moving to conduct a BINNED likelihood analysis."

    obs = BinnedObs(binnedExpMap= Name + '_expcube_lc.fits',expCube= Name + '_ltcube_lc.fits',srcMaps= Name + '_srcmaps_lc.fits', irfs='P8R2_SOURCE_V6')
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
        limit,results = IUL.calc_int(analysis, Name,cl=0.95)

    #Run gtselect to make smaller data fits file to compute the exposure, set to 3 degrees around source of interest
    f._setRad(3)
    f._setInfile(Evfile)
    f._setOutfile( Name + '_gtselect_exposure.fits')
    print "Creating file " + Name + "_gtselect_exposure.fits"
    f.amonSelect()

    #Run gtmaketime on this small region
    f._setEvfile( Name + '_gtselect_exposure.fits')
    f._setOutfile( Name + '_gtmktime_exposure.fits')
    print('Working on file ' + str(f.getOutfile()))
    f.amonTime()

    my_apps.evtbin['algorithm'] = 'LC'
    my_apps.evtbin['evfile'] = f.getOutfile()
    my_apps.evtbin['outfile'] = Name + '_LC_exposure.fits'
    my_apps.evtbin['scfile'] = f.getScfile()
    my_apps.evtbin['tbinalg'] = 'LIN'
    my_apps.evtbin['tstart'] = f.getTmin()
    my_apps.evtbin['tstop'] = f.getTmax()
    my_apps.evtbin['dtime'] = TSTOP - TSTART
    my_apps.evtbin.run()

    yes = subprocess.call(['gtexposure', Name + '_LC_exposure.fits',f.getScfile(),'P8R2_SOURCE_V6',Name + '_output_model.xml',Name])
    if yes == 0:
        print "Exposure map has been created"
    else:
        print "Subprocessing failed. Unable to create exposure map with gtexposure."
    
    print "Time bin complete."

    hdulist = pyfits.open(Name + '_LC_exposure.fits')
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

"""

		This portion of the code will split the data into bands and run a likelihood analysis
		from start to finish on the data.

"""
def run(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,TSTOP,Evfile,bins,zmax,evclass,evtype,TSul,NMtol, bin_num,runMRM):
	TTime = TSTOP - TSTART
	bin_size = TTime/bin_num
	i = 0
	while i < bin_num:
		tstop = TSTART + bin_size
		runFermiTools(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,tstop,Evfile,bins,zmax,evclass,evtype,TSul,NMtol, i + 1,runMRM)
		i = i + 1
		TSTART = TSTART + bin_size
	print 'LC data has been compiled. Please see lc_output.txt for results.'
"""
    Copy of makeResMap - needed to copy to avoid error from farith everytime residuals map was made
"""
def runMakeResMap(Name, Ra, Dec, radius, binsz,minEnergy,maxEnergy,lc_bin_num):
    print "This is makeResMap.\nStarting construction of residuals map.\n"
    f = FermiObject()
    print "Creating model map of the region surrounding " + str(Name) + " with gtmodel."
    f._setSrcmaps(Name + '_srcmaps_lc.fits')
    f._setExpcube(Name + '_ltcube_lc.fits')
    f._setBexpmap(Name + '_expcube_lc.fits')
    f._setSrcmdl(Name + '_output_model.xml')
    f._setOutfile( Name + '_model_map_' + str(lc_bin_num) + '.fits')
    f._setIrfs('CALDB')
    f.amonModelmaps()
    print "Finished model map. Now constructing counts map for comparison."
    f._setAlgorithm('CMAP')
    f._setEvfile(Name + '_gtmktime_lc.fits')
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

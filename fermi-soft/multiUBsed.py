"""

        multiUBsed.py - creates a spectrum for a source of interest
        using an unbinned likelihood analysis. Be sure to have a spectral model for
        your source of interest in your output model from your likelihood analysis
        for the entire dataset; you will be using this XML file throughout.

        Author: Michael Toomey <mwt5345@psu.edu>
        Institution: The Pennsylvania State University

        * Relies on the implementation of the FermiObject.py to gain access to the
            use of the various Fermi Tools required for this analysis.

"""

from FermiObject import FermiObject
import sys
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi/gtapps_mp-1.1')
from UnbinnedAnalysis import *
import IntegralUpperLimit as IUL
import numpy as np
import subprocess
from math import sqrt,log10
import gt_apps as my_apps
import pyfits
import multiprocessing
import pyLikelihood

# Conducts the complete likelihood analysis to get the flux, flux error, and
# TS value for each bin of the light curve as well as upper limits when needed.
def runFermiTools(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,TSTOP,Evfile,bins,zmax,evclass,evtype,TSul,NMtol,sed_bin_num,runMRM):

    if evclass == 512:
        irf = "P8R2_ULTRACLEAN_V6"
    elif evclass == 128:
        irf = "P8R2_SOURCE_V6"
    elif evclass == 256:
        irf = "P8R2_CLEAN_V6"
    elif evclass == 1024:
        irf = "P8R2_ULTRACLEANVETO_V6"

    print "Working on bin " + str(sed_bin_num) + " for the spectrum."

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
    f._setOutfile( Name + '_gtselect' + str(sed_bin_num) + '_sed.fits')
    f.amonSelect()
    print('File cuts have been made. Now making cuts for GTI using spacecraft file.')

    """

        Following steps execute Fermi Tool gtmktime

    """

    f._setScfile(SCFile)
    f._setRoicut('no')
    f._setEvfile( Name + '_gtselect' + str(sed_bin_num) + '_sed.fits')
    f._setOutfile( Name + '_gtmktime' + str(sed_bin_num) + '_sed.fits')
    ###############################################
    #         Filter expression                   #
    Filter = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
    ###############################################
    f._setFilter(Filter)
    print('Working on file ' + str(f.getOutfile()) + '. . .')
    f.amonTime()
    print('File cuts have been made.')
    print('Using XML model from whole dataset.\n Moving on to gtltcube.')

    """ Livetime need only be calculated once, not dependant on model, only dependant on Fermi orientation
    print "Now working on ltcube file using gtltcube\n"
    my_apps.expCube['evfile'] = Name + '_gtmktime' + str(sed_bin_num) + '_sed.fits'
    my_apps.expCube['scfile'] =  SCFile
    my_apps.expCube['outfile'] = Name + '_ltcube' + str(sed_bin_num) + '_sed.fits'
    my_apps.expCube['dcostheta'] = 0.025
    my_apps.expCube['binsz'] = 1
    my_apps.expCube['phibins'] = 0
    my_apps.expCube['zmax'] = zmax
    my_apps.expCube['chatter'] = 0
    my_apps.expCube.run()"""

    print "\nltcube complete.\nMoving to compute exposure map with gtexpmap.\n"
    my_apps.expMap['evfile'] = Name + '_gtmktime' + str(sed_bin_num) + '_sed.fits'
    my_apps.expMap['scfile'] = SCFile
    my_apps.expMap['expcube'] = Name + '_ltcube.fits'# + str(sed_bin_num) + '_sed.fits'
    my_apps.expMap['outfile'] = Name + '_expMap' + str(sed_bin_num) + '_sed.fits'
    my_apps.expMap['irfs'] = 'CALDB'
    my_apps.expMap['srcrad'] = radius + 10
    my_apps.expMap['nlong'] = 4*(radius + 10)
    my_apps.expMap['nlat'] = 4*(radius + 10)
    ebin = int(10*log10(maxEnergy/minEnergy))
    print "There are " + str(ebin) + " energy bands."
    my_apps.expMap['nenergies'] = ebin
    my_apps.expMap.run()
    print "Finnished making exposure map.\n"

    print "Calcualting the diffuse response for photons in this bin."    
    my_apps.diffResps['evfile'] = Name + '_gtmktime' + str(sed_bin_num) + '_sed.fits'
    my_apps.diffResps['scfile'] = SCFile
    my_apps.diffResps['srcmdl'] = Name + '_output_model.xml'
    my_apps.diffResps['irfs'] = 'CALDB'
    my_apps.diffResps.run()

    print "Finished calculating diffuse response. Now moving to conduct a BINNED likelihood analysis."

    obs = UnbinnedObs(Name + '_gtmktime' + str(sed_bin_num) + '_sed.fits', SCFile,expMap= Name + '_expMap' + str(sed_bin_num) + '_sed.fits',expCube= Name + '_ltcube.fits', irfs=irf)
    analysis = UnbinnedAnalysis(obs,Name + '_output_model.xml', optimizer='NewMinuit')
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
    ESUM = minEnergy + maxEnergy
    EMID = ESUM/2
    limit = Flux
    if analysis.Ts(Name) < TSul:
        UL = True
        limit,results = IUL.calc_int(analysis, Name,cl=0.90, emin=minEnergy, emax=maxEnergy)
    """
    #Run gtselect to make smaller data fits file to compute the exposure, set to 3 degrees around source of interest
    f._setRad(3)
    f._setInfile(Evfile)
    f._setOutfile( Name + '_gtselect' + str(sed_bin_num) + '_exposure.fits')
    print "Creating file " + Name + "_gtselect_exposure.fits"
    f.amonSelect()

    #Run gtmaketime on this small region
    f._setEvfile( Name + '_gtselect' + str(sed_bin_num) + '_exposure.fits')
    f._setOutfile( Name + '_gtmktime' + str(sed_bin_num) + '_exposure.fits')
    print('Working on file ' + str(f.getOutfile()))
    f.amonTime()

    my_apps.evtbin['algorithm'] = 'LC'
    my_apps.evtbin['evfile'] = f.getOutfile()
    my_apps.evtbin['outfile'] = Name + '_LC' + str(sed_bin_num) + '_exposure.fits'
    my_apps.evtbin['scfile'] = f.getScfile()
    my_apps.evtbin['tbinalg'] = 'LIN'
    my_apps.evtbin['tstart'] = f.getTmin()
    my_apps.evtbin['tstop'] = f.getTmax()
    my_apps.evtbin['dtime'] = TSTOP - TSTART
    my_apps.evtbin.run()

    yes = subprocess.call(['gtexposure', Name + '_LC' + str(sed_bin_num) + '_exposure.fits',f.getScfile(),irf,Name + '_output_model.xml',Name])
    if yes == 0:
        print "Exposure map has been created"
    else:
        print "Subprocessing failed. Unable to create exposure map with gtexposure."
    
    print "Time bin complete."

    hdulist = pyfits.open(Name + '_LC' + str(sed_bin_num) + '_exposure.fits')
    tbdata = hdulist[1].data
    z = tbdata['EXPOSURE']
    exp = z[0]
    """
    ################################################################
    #           This portion prints to the text file               #
    ################################################################
    f = open("sed_output.txt", "a")
    f.write(str(Flux) + ',' + str(Ferr) + ',' + str(ef) + ',' + str(ef_err) + ',' + str(limit) + ','+ str(analysis.Ts(Name)) + ',' + str(UL) + ',' + str(EMID) + '\n')
    f.close()
    print "Likelihood analysis on this band is complete."

    yes = subprocess.call(['rm',Name + '_gtselect' + str(sed_bin_num) + '_sed.fits',Name + '_gtmktime' + str(sed_bin_num) + '_sed.fits',Name + '_cmap' + str(sed_bin_num) + '_sed.fits',Name + '_ccube' + str(sed_bin_num) + '_sed.fits',Name + '_expMap' + str(sed_bin_num) + '_sed.fits',Name + '_srcmaps' + str(sed_bin_num) + '_sed.fits'])
    if yes == 0:
        print 'Files for bin have been deleted'
    else:
        print "Subprocessing failed. Unable to delete files for bin."

"""

        This portion of the code will split the data into bands and run a likelihood analysis
        from start to finish on the data.

"""
def run(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,binsz,TSTART,TSTOP,Evfile,bins,zmax,evclass,evtype,TSul,NMtol, bin_num,runMRM,num):
    f = open("sed_output.txt", "a")
    f.write("#Flux,FluxErr,EnFLux,EnFluxErr,upperlimit?,TS,logMidEn\n")
    f.close()
    x0 = log10(minEnergy)
    x1 = log10(maxEnergy)
    dx = x1 - x0
    delx = dx/num
    i = 0
    pool = multiprocessing.Pool(num)
    while i < bin_num:
        emax = 10 ** (x0 + delx)
        emin = 10 ** x0 
        pool.apply_async(runFermiTools, args=(Name,RA,DEC,emin,emax,SCFile,radius,binsz,TSTART,TSTOP,Evfile,bins,zmax,evclass,evtype,TSul,NMtol, i + 1,runMRM))
        i = i + 1
        x0 = x0 + delx
    pool.close()
    pool.join()
    print 'LC data has been compiled. Please see lc_output.txt for results.'

"""

	multiLike.py - this script picks up from where makeFermiFiles.py
	left off, if you are planing on conducting an BINNNED Likelihood analysis.
	Will produce upper limits. This is similar to fromXMLtolike.py with
	the exception that it runs two fits on the data, one loose fit
	is initially done with MINUIT then a much tighter fit is done with NEWMINUIT.
    You also have the option to run a single fit with NEWMINUIT by using the fucntion
    miniLike. Both of the functions have the capability to make a residual map. Both
    multiLike and miniLike functions have the capability to spawn two jobs when both
    a residuals map and upper limit are needed.

	Author: Michael Toomey <mwt5345@psu.edu>
	Institution: The Pennsylvania State University

	This version has been specifically made for use on the NRL's HESELIN system.


"""
import sys
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi/gtapps_mp-1.1')
sys.path.insert(0,'/storage/home/mwt5345/work/python/Fermi')
from FermiObject import FermiObject
from BinnedAnalysis import *
import numpy as np
import IntegralUpperLimit as IUL    #Bayesian
from UpperLimits import UpperLimits #Frequentist
import subprocess
from math import sqrt
from makeResMap import runMakeResMap
from multiprocessing import Process
from ConfigParser import ConfigParser

"""

	Runs Fermi Tool gtsrcmaps

"""
def runAmonSrcmaps(Name,SCFile):
	
	f = FermiObject()

	print "Now begining computation of model counts map with gtsrcmaps.\n"

	f._setScfile(SCFile)
	f._setExpcube( '../' + Name + '_ltcube.fits')
	f._setCmap( '../' + Name + '_ccube.fits')
	f._setBexpmap( '../' + Name + '_expcube.fits')
	f._setOutfile( Name + '_srcmaps.fits')
	f._setIrfs('CALDB')
	f._setSrcmdl(Name + '_in_model.xml')
	f.amonSrcmaps()
	
	print "Finished srcmaps. Now moving to conduct a BINNED likelihood analysis."

"""

	Runs a low tolerance likelihood analysis with the optimizer MINUIT

"""
def runMINUIT(obs, Name, Mtol):
	like1 = BinnedAnalysis(obs,srcModel= Name + '_in_model.xml', optimizer='MINUIT')
	like1.tol = Mtol
	like1obj = pyLike.Minuit(like1.logLike)
	like1.fit(verbosity=3,covar=True,optObject=like1obj)
	like1.saveCurrentFit()
	print "Minuit has converged whith Code " + str(like1obj.getQuality())
	like1.logLike.writeXml( Name + '_fit1.xml')
	return like1

"""

	Runs a low tolerance likelihood analysis with the optimizer NEWMINUIT

"""
def runNEWMINUITonly(obs, Name,NMtol):
	print "\nNow begining fit with NEWMINUIT"
	analysis = BinnedAnalysis(obs,srcModel= Name + '_in_model.xml', optimizer='NEWMINUIT')
	likeObj = pyLike.NewMinuit(analysis.logLike)
	analysis.tol = NMtol
	analysis.fit(verbosity=3,covar=True,optObject=likeObj)
	analysis.saveCurrentFit()
	print "Likelihood has converged whith Code " + str(likeObj.getRetCode())
	analysis.writeXml( Name + '_output_model.xml')
	return analysis

def runNEWMINUIT(obs, Name,NMtol):
	print "\nNow begining new tighter fit with NEWMINUIT"
	analysis = BinnedAnalysis(obs,srcModel= Name + '_fit1.xml', optimizer='NEWMINUIT')
	likeObj = pyLike.NewMinuit(analysis.logLike)
	analysis.tol = NMtol
	analysis.fit(verbosity=3,covar=True,optObject=likeObj)
	analysis.saveCurrentFit()
	print "Likelihood has converged whith Code " + str(likeObj.getRetCode())
	analysis.writeXml( Name + '_output_model.xml')
	return analysis

"""
	
	Prints the results of the likelihood analysis

"""
def printResults(analysis, Name, minEn, maxEn):
    MeVtoErg = 1.602e-6
    print "\nNow displaying results of the fit.\n"
    print "+---------------------------------------------------------------------------------+"
    for source in analysis.sourceNames():
        if(np.shape(analysis.freePars(source))[0] > 0):
            print "******This was a free source.*******"
            print "+---------------------------------------------------------------------------------+"
            print analysis.model[source]
            print "TS value: " +  str(analysis.Ts(source))
            print "Flux: " + str(analysis.flux(source, emin=minEn, emax=maxEn)) + " +/- " + str(analysis.fluxError(source, emin=minEn, emax=maxEn)) + " photons/cm^2/s"
            print "0.1 - 100 GeV Flux: ", analysis.flux(source, emin=100, emax=1e5), " +/- ", analysis.fluxError(source, emin=100, emax=1e5), " photons/cm^2/s"
            print "Energy Flux: " + str(analysis.energyFlux(source,minEn,maxEn)*MeVtoErg) + " +/- " + str(analysis.energyFluxError(source,minEn,maxEn)*MeVtoErg) + " erg/cm^2/s"
            print "0.1 - 100 GeV Energy Flux: " + str(analysis.energyFlux(source,1e2,1e5)*MeVtoErg) + " +/- " + str(analysis.energyFluxError(source,1e2,1e5)*MeVtoErg) + " erg/cm^2/s"
            print "Predicted number of counts: " + str(analysis.NpredValue(source))
            print "+---------------------------------------------------------------------------------+"
        else:
            pass
    print "+---------------------------------------------------------------------------------+"
    print "\nResults for source of interest: " + str(analysis.model[Name])
    print "\nFlux: ", analysis.flux(Name, emin=minEn, emax=maxEn), " +/- ", analysis.fluxError(Name, emin=minEn, emax=maxEn), " photons/cm^2/s"
    print "0.1 - 100 GeV Flux: ", analysis.flux(Name, emin=100, emax=1e5), " +/- ", analysis.fluxError(Name, emin=100, emax=1e5), " photons/cm^2/s"
    print "Energy Flux: " + str(analysis.energyFlux(Name,minEn,maxEn)*MeVtoErg) + " +/- " + str(analysis.energyFluxError(Name,minEn,maxEn)*MeVtoErg) + " erg/cm^2/s"
    print "0.1 - 100 GeV Energy Flux: " + str(analysis.energyFlux(Name,1e2,1e5)*MeVtoErg) + " +/- " + str(analysis.energyFluxError(Name,1e2,1e5)*MeVtoErg) + " erg/cm^2/s"
    print "TS Value: ", analysis.Ts(Name,reoptimize=False)
    print "Apx. SD: ", np.sqrt(analysis.Ts(Name))
    print "\nLikelihood analysis is complete."
    print "+---------------------------------------------------------------------------------+"
    print "+---------------------------------------------------------------------------------+"
    #Writes out results.dat and counts_spectra.fits file with results of the fit for all sources.
    analysis.writeCountsSpectra()
"""

	Calculates the Bayesian upper limit of a source to 95% confidence

"""
def runBayesianUL(analysis, Name, minEnergy, maxEnergy):
    print "\nThe TS is below the threshold, calculating 95% confidence-level Bayesian upper limit."
    limit,results = IUL.calc_int(analysis,Name,cl=0.95,emin=minEnergy, emax=maxEnergy)
    print "Bayesian upper limit: " + str(limit) + " photons/cm^2/s"
    if minEnergy != 100 and maxEnergy != 1e5:
        print "Calculating greater than 100 MeV upper limit flux"
        limit,results = IUL.calc_int(analysis,Name,cl=0.95,emin=100,emax=1e5)
        print "0.1 - 100 GeV Bayesian upper limit: " + str(limit) + " photons/cm^2/s"

"""

	Calculates an upperlimit using the frequentist approach

"""
def runFreqUL(analysis, Name, minEnergy, maxEnergy):
    print "\nThe TS is below the threshold, calculating 95% confidence-level upper limit with frequentist approach."
    ul = UpperLimits(analysis)
    #The delta argument is used to determin the confidence level of upper limit. Derived from the chi distribution, 3.841 corresponds
    #to a p-value of 0.05 for 1 degree of freedom (here the fit of the normilizations) which can be interpreted as a 95% upper limit.
    #The value is divided by two because 2*logLikelihood difference = chi distribution value.
    ul[Name].compute(emin=minEnergy, emax=maxEnergy, delta=3.841/2.)
    print "Results:"
    print ul[Name].results
    if minEnergy != 100 and maxEnergy != 1e5:
        print "Calculating 0.1 - 100 GeV upper limit flux"
        ul[Name].compute(emin=100, emax=1e5, delta=3.841/2.)
        print "Results:"
        print ul[Name].results

def checkValues(ulTS, BUL, SRC):
	if isinstance(ulTS, (int,float)):
		pass
	else:
		raise ValueError('Invalid TS value. Please enter a float or int.')
	if isinstance(BUL, bool):
		pass
	else:
		raise ValueError('Invalid value. Please enter bool True for Bayesian UL or False for Freq. UL.')
	if isinstance(SRC, bool):
		pass
	else:
		raise ValueError('Invalid value. Please enter bool True to compute source map, False to skip.')

"""

	This fuction computes the likelihood, you have the ability to either calculate
	the source maps or skip the calculation and go straight to the likelihood analysis.
	The user also has the ability to set upper limits starting at a TS of their choice
	and using either a Bayesian or frequentist approach. Will produce a residuals map
    if prompted.

	Name - name of your source as a string
	Ra - right ascension of center of ROI as int/float
	Dec - declination of center of ROI as int/float
	SCFile - the path to the directory of your spacecraft file as String
	ulTS - value of TS below which your source of interest will have upper limit calculated as int/float
	BUL - True for Bayesian upper limits, Fasle for Frequentist as bool
	SRC - True to make source maps, False to skip as bool
    Mtol - tolerance for MINUIT as int/float
    NMtol - tolerance for NEWMINUIT as int/float
    radius - the radius of your region of interest in degrees as int/float
    binsz - the size of each pixels in degrees as int/float
    makeRes - True to make residual map, False to skip as bool

"""
def multiLike():

    parser = ConfigParser()
    parser.read('../config.ini')

    #parameters from config file
    Name = str(parser.get('params','Name'))
    Ra = float(parser.get('params','RA'))
    Dec = float(parser.get('params','DEC'))
    minEnergy = float(parser.get('params','minEnergy'))
    maxEnergy = float(parser.get('params','maxEnergy'))
    SCFile = str(parser.get('params','SCFile'))
    radius = float(parser.get('params','radius'))
    binsz = float(parser.get('params','binsz'))
    ulTS = float(parser.get('params','ulTS'))
    BUL = bool(parser.get('params','BUL'))
    SRC = bool(parser.get('params','SRC'))
    Mtol = float(parser.get('params','Mtol'))
    NMtol = float(parser.get('params','NMtol'))
    makeRes = bool(parser.get('params','makeRes'))
    evclass = int(parser.get('params','evclass'))

    if evclass == 512:
        irf = "P8R2_ULTRACLEAN_V6"
    elif evclass == 128:
        irf = "P8R2_SOURCE_V6"
    elif evclass == 256:
        irf = "P8R2_CLEAN_V6"
    elif evclass == 1024:
        irf = "P8R2_ULTRACLEANVETO_V6"
    print "This is multiLike.\nPlease make sure that you have added a model for your source with the name: " + str(Name)
	#Checks that values not used until later in the program have the correct type
    checkValues(ulTS, BUL, SRC)

	#Calls function to compute the source map with using the input .xml model
    if SRC == True:
        runAmonSrcmaps(Name,SCFile)

	#Likelihood analysis is done with a looser tolerance with MINUIT first, followed by a tight tolerance with NEWMINUIT
    obs = BinnedObs(binnedExpMap= '../' + Name + '_expcube.fits',expCube= '../' + Name + '_ltcube.fits',srcMaps= Name + '_srcmaps.fits', irfs=irf)
	
	#The first fit is done with Minuit 
    runMINUIT(obs,Name,Mtol)

	#The second fit is done with NEWMINUIT
    analysis = runNEWMINUIT(obs,Name,NMtol)

	#Will print out the TS and flux values for all the sources
    printResults(analysis, Name, minEnergy, maxEnergy)

    ##########################################################################
    ###                      Multiprocessing bit                           ###
    ##########################################################################
    if  (makeRes == True) and (BUL == True) and (analysis.Ts(Name) < ulTS):
        print "Breaking into two threads . . . . . . . . ."
        p1 = Process(target=runMakeResMap,args=(Name, Ra, Dec, minEnergy,maxEnergy,radius, binsz))
        print "Starting residuals map thread . . . . . . ."
        p1.start()
        p2 = Process(target=runBayesianUL,args=(analysis,Name,minEnergy,maxEnergy))
        print "Starting upper limit thread . . . . . . . ."
        p2.start()
        p1.join()
        p2.join()
    elif (makeRes == True) and (BUL == False) and (analysis.Ts(Name) < ulTS):
        print "Breaking into two threads . . . . . . . . ."
        p1 = Process(target=runMakeResMap,args=(Name,Ra,Dec,minEnergy,maxEnergy,radius,binsz))
        print "Starting residuals map thread . . . . . . ."
        p1.start()
        p2 = Process(target=runFreqUL,args=(analysis,Name,minEnergy,maxEnergy))
        print "Starting upper limit thread . . . . . . . ."
        p2.start()
        p1.join()
        p2.join()
    elif  (makeRes == False) and (analysis.Ts(Name) < ulTS):
        if BUL == True:
            runBayesianUL(analysis,Name,minEnergy,maxEnergy)
        else:
            runFreqUL(analysis,Name,minEnergy,maxEnergy)
    else:
        if makeRes == True:
            runMakeResMap(Name, Ra, Dec, minEnergy,maxEnergy,radius, binsz)
        else:
            pass

    print "\nProgram Finished."

def miniLike():

    parser = ConfigParser()
    parser.read('../config.ini')

    #parameters from config file
    Name = str(parser.get('params','Name'))
    Ra = float(parser.get('params','RA'))
    Dec = float(parser.get('params','DEC'))
    minEnergy = float(parser.get('params','minEnergy'))
    maxEnergy = float(parser.get('params','maxEnergy'))
    SCFile = str(parser.get('params','SCFile'))
    radius = float(parser.get('params','radius'))
    binsz = float(parser.get('params','binsz'))
    ulTS = float(parser.get('params','ulTS'))
    BUL = bool(parser.get('params','BUL'))
    SRC = bool(parser.get('params','SRC'))
    NMtol = float(parser.get('params','NMtol'))
    makeRes = bool(parser.get('params','makeRes'))
    evclass = int(parser.get('params','evclass'))

    if evclass == 512:
        irf = "P8R2_ULTRACLEAN_V6"
    elif evclass == 128:
        irf = "P8R2_SOURCE_V6"
    elif evclass == 256:
        irf = "P8R2_CLEAN_V6"
    elif evclass == 1024:
        irf = "P8R2_ULTRACLEANVETO_V6"
    print "This is multiLike.\nPlease make sure that you have added a model for your source with the name: " + str(Name)
	#Checks that values not used until later in the program have the correct type
    checkValues(ulTS, BUL, SRC)

	#Calls function to compute the source map with using the input .xml model
    if SRC == True:
        runAmonSrcmaps(Name, SCFile)
    else:
        pass

	#Likelihood analysis is done with a looser tolerance with MINUIT first, followed by a tight tolerance with NEWMINUIT
    obs = BinnedObs(binnedExpMap= '../' + Name + '_expcube.fits',expCube= '../' + Name + '_ltcube.fits',srcMaps= Name + '_srcmaps.fits', irfs=irf)
	
	#The first fit is done with Minuit 
    analysis = runNEWMINUITonly(obs,Name,NMtol)
	
	#Will print out the TS and flux values for all the sources
    printResults(analysis, Name, minEnergy, maxEnergy)

    ##########################################################################
    ###                      Multiprocessing bit                           ###
    ##########################################################################
    if  (makeRes == True) and (BUL == True) and (analysis.Ts(Name) < ulTS):
        print "Breaking into two threads . . . . . . . . ."
        p1 = Process(target=runMakeResMap,args=(Name, Ra, Dec, minEnergy,maxEnergy,radius, binsz))
        print "Starting residuals map thread . . . . . . ."
        p1.start()
        p2 = Process(target=runBayesianUL,args=(analysis,Name,minEnergy,maxEnergy))
        print "Starting upper limit thread . . . . . . . ."
        p2.start()
        p1.join()
        p2.join()
    elif (makeRes == True) and (BUL == False) and (analysis.Ts(Name) < ulTS):
        print "Breaking into two threads . . . . . . . . ."
        p1 = Process(target=runMakeResMap,args=(Name,Ra,Dec,minEnergy,maxEnergy,radius,binsz))
        print "Starting residuals map thread . . . . . . ."
        p1.start()
        p2 = Process(target=runFreqUL,args=(analysis,Name,minEnergy,maxEnergy))
        print "Starting upper limit thread . . . . . . . ."
        p2.start()
        p1.join()
        p2.join()
    elif  (makeRes == False) and (analysis.Ts(Name) < ulTS):
        if BUL == True:
            runBayesianUL(analysis,Name,minEnergy,maxEnergy)
        else:
            runFreqUL(analysis,Name,minEnergy,maxEnergy)
    else:
        if makeRes == True:
            runMakeResMap(Name, Ra, Dec, minEnergy,maxEnergy,radius, binsz)
        else:
            pass

    print "\nProgram Finished."

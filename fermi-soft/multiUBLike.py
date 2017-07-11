"""

	multiUBLike.py - this script picks up from where makeUBFermiFiles.py
	left off, if you are planing on conducting an UNBINNNED Likelihood analysis.
	Will produce upper limits.

	Author: Michael Toomey <mwt5345@psu.edu>
	Institution: The Pennsylvania State University

"""
import multiLike
import gt_apps as my_apps
from UnbinnedAnalysis import *
import IntegralUpperLimit as IUL
import numpy as np
import subprocess
from math import sqrt,log10
import gt_apps as my_apps
import pyfits
import multiprocessing
import pyLikelihood

def run(Name, Ra, Dec, minEnergy, maxEnergy, SCFile, ulTS, NMtol,evclass,model):

    if evclass == 512:
        irf = "P8R2_ULTRACLEAN_V6"
    elif evclass == 128:
        irf = "P8R2_SOURCE_V6"
    elif evclass == 256:
        irf = "P8R2_CLEAN_V6"
    elif evclass == 1024:
        irf = "P8R2_ULTRACLEANVETO_V6"

    print "This is multiUBLike.\nPlease make sure that you have added a model for your source with the name: " + str(Name)

    print "Calcualting the diffuse response for photons in this bin."    
    my_apps.diffResps['evfile'] = '../' + Name + '_gtmktime.fits'
    my_apps.diffResps['scfile'] = SCFile
    my_apps.diffResps['srcmdl'] = model
    my_apps.diffResps['irfs'] = irf
    my_apps.diffResps.run()

    print "Finished calculating diffuse response. Now moving to conduct a UNBINNED likelihood analysis."

    obs = UnbinnedObs('../' + Name + '_gtmktime.fits', SCFile,expMap= '../' + Name + '_expMap.fits',expCube= '../' + Name + '_ltcube.fits', irfs=irf)
    analysis = UnbinnedAnalysis(obs,model,optimizer='NewMinuit')
    likeObj = pyLike.NewMinuit(analysis.logLike)
    analysis.tol = NMtol
    lkl = analysis.fit(verbosity=0,covar=True,optObject=likeObj)
    analysis.writeXml( Name + '_output_model.xml')
    fit = likeObj.getRetCode()
    print "Likelihood has converged whith Code " + str(likeObj.getRetCode())

    multiLike.printResults(analysis,Name, minEnergy,maxEnergy)
    print "Fit has likelihood: " + str(lkl)

    print "\nThe TS is below the threshold, calculating 95% confidence-level Bayesian upper limit."
    limit,results = IUL.calc_int(analysis,Name,cl=0.95,emin=minEnergy, emax=maxEnergy)
    print "Bayesian upper limit: " + str(limit) + " photons/cm^2/s"
    #Array that returns the results of the unbinned analysis
    #log-likelihood,flux,flux_err,test statisitc
    Return = [lkl,analysis.flux(Name, emin=minEn, emax=maxEn),analysis.fluxError(Name, emin=minEn, emax=maxEn),limit,analysis.Ts(Name,reoptimize=False)]
    return Return

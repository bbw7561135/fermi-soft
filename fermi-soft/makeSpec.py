from bdlikeSED_v14 import *
from BinnedAnalysis import *
import sys
sys.path.insert(0,'/home/tyrelj/my_python')

"""
    makeSpec.py - implements tyrels code bdlikeSED_14 to make a spectral plot

    Name    # The name of your source <string>
    Bin     # The number of bins for your SED plot <int>
    RadLim  # The radius out to which sources above SigLim will be set free <int/float>
    tslim   # The TS limit for plotting a point vs. plotting an upper limit <int>
    nplim   # The predicted number of counts limit for plotting a point vs. plotting an upper limit <int>
    plot    # True to plot, false to skip <bool>

    # For low TS source set tslim=-1 and nplim=0

"""
def run(Name,Bin,RadLim,SigLim,tslim,nplim,plot):
    obs=BinnedObs(srcMaps= Name + '_srcmaps.fits', expCube= '../' + Name + '_ltcube.fits', binnedExpMap= '../' + Name + '_expcube.fits', irfs='CALDB')
    like=BinnedAnalysis(obs, Name + '_output_model.xml', 'NEWMINUIT')
    inputs=bdlikeInput(like, Name, '../' + Name + '_gtmktime.fits', nbins=Bin) 
    inputs.plotBins()
    inputs.getfullFit(radLim=RadLim,sigLim=SigLim)
    sed=bdlikeSED(inputs)
    sed.logMidECent()
    sed.fitBands()
    sed.Plot(tslim,nplim,plot)

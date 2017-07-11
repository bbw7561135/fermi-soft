import sys
sys.path.insert(0,'/Projects/GLAST/share/bin')
sys.path.insert(0,'/home/mtoomey/FilesAndScripts')
from math import log10,sqrt
import gt_apps as my_apps
import subprocess

def run(Name,RA,DEC,minEnergy,maxEnergy,SCFile,radius,model):
    print "This is findSource v1.0"
    #gtexpmap
    my_apps.expMap['evfile'] = '../' + Name + '_gtmktime.fits'
    my_apps.expMap['scfile'] = SCFile
    my_apps.expMap['expcube'] = '../' + Name + '_ltcube.fits'
    my_apps.expMap['outfile'] = Name + '_expmap_ub.fits'
    my_apps.expMap['irfs'] = 'CALDB'
    my_apps.expMap['srcrad'] = radius + 10
    my_apps.expMap['nlong'] = 4*(radius + 10)
    my_apps.expMap['nlat'] = 4*(radius + 10)
    ebin = int(10*log10(maxEnergy/minEnergy))
    print "There are " + str(ebin) + " energy bans."
    my_apps.expMap['nenergies'] = ebin
    my_apps.expMap.run()

    yes = subprocess.call(['gtfindsrc evfile=../' + str(Name) + '_gtmktime.fits scfile=' + str(SCFile) + ' outfile=' + str(Name) + '_position irfs=CALDB expcube=../' + str(Name) + '_ltcube.fits expmap=' + str(Name) + '_expmap_ub.fits srcmdl=' + str(model) + ' target=' + str(Name) + ' coordsys=CEL ra=' + str(RA) + ' dec=' + str(DEC) + ' optimizer=MINUIT ftol=1e-2 atol=0.01'],shell=True)
    if yes == 0:
        print "Best fit source position has been found, see " + str(Name) + "_position"
    else:
        print "Subprocessing failed. Unable to find source position with gtfindsource."

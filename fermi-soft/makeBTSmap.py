"""

	makeBTSmap.py - makes a binned TS map using Fermi Tools
	python capabilities.

	Author: Michael Toomey <mwt5345@psu.edu>
	Institution: The Pennsylvania State University

	This version was specifically made to operate on the 
	HESELIN system at the NRL in conjuction with the
	scripts createXMLfromscratch.py and fromXMLtolike.py

"""

import gt_apps as my_apps
import pyfits

"""

    Name =      # The name of your source here <String>
    Ra =        # The right ascension of your source here <int/float>
    Dec =       # The declination of your source here <int/float>
    SCFILE =    # The complete path to your spacecraft file <String>
    model =     # The name of you output source model to use <String>
    pix =       # The number of pixels on a side i.e. make a TS map of pix by pix size <int/float>
    pixsz =     # The size of your pixels <int/float>
    ftol =      # The tolerance of the likelihood fit <int/float>

"""
def run(Name, Ra, Dec, SCFILE,model,pix,pixsz,ftol):
    #Does the equivelent of the Fermi Tool gttsmap
    my_apps.TsMap['statistic'] = "BINNED"
    my_apps.TsMap['scfile'] = "/Projects/GLAST/data/FT2-dev.fits"
    my_apps.TsMap['cmap'] = "../" + Name + "_ccube.fits"
    my_apps.TsMap['evfile'] = "../" + Name + "_gtmktime.fits"
    my_apps.TsMap['bexpmap'] = "../" + Name + "_expcube.fits"
    my_apps.TsMap['expcube'] = "../" + Name + "_ltcube.fits"
    my_apps.TsMap['srcmdl'] = model
    my_apps.TsMap['irfs'] = "CALDB"
    my_apps.TsMap['optimizer'] = "NEWMINUIT"
    my_apps.TsMap['outfile'] = Name + "_tsmap.fits"
    my_apps.TsMap['nxpix'] = pix
    my_apps.TsMap['nypix'] = pix
    my_apps.TsMap['binsz'] = pixsz
    my_apps.TsMap['coordsys'] = "CEL"
    my_apps.TsMap['xref'] = Ra
    my_apps.TsMap['yref'] = Dec
    my_apps.TsMap['proj'] = 'AIT'
    my_apps.TsMap['ftol'] = ftol
    my_apps.TsMap.run()

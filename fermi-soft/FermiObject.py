"""
    FermiObject.py - python class that creates a python object for handling Fermi
        LAT data analysis. Utilizes Fermi Tools python scripts.
        
    Author: Michael Toomey <mwt5345@psu.edu>
    Institution: The Pennsylvania State University
    Group: The Astrophysical Multimessenger Observatory Network
    Last Updated: June 21, 2016 10:46:34 EDT

    *Due to a recent increase in capability of this class, not all
     safeguards have been put in place yet.
"""

import gt_apps as my_apps

class FermiObject(object):

    """
    Instances of class including:
            sourcename, evclass,evtype,ra,dec,rad,emin,emax,zmax,tmin,tmax,infile,
         outfile, evfile, scfile, expcube, dcostheta, binsz, irfs, srcrad, nlong,
         nlat, nenergies, srcmdl, roicut and many more.
        Each instance has a corresponding get<name>() and _set<name> module to access 
     and change objects instance values.
        
        Example: 
            >>> fermi1 = FermiObject()
            >>> fermi1.getEvclass()
            AttributeError: 'trimData' object has no attribute '_evclass'
            >>> fermi1._setEvclass(128)
            >>> fermi1.getEvclass()
            128
    """

    print "Fermi Science Tools are being implemented via FermiObject\nThis version of FermiObject was made for Science Tools 10-00-05"

    def getSourcename(self):
        return self._sourcename

    def _setSourcename(self, sn):
          if isinstance(sn, basestring):
               self._sourcename = sn
          else:
               raise ValueError('Invalid source name. Please enter source name or event number.')
        
    def getEvclass(self):
        return self._evclass
        
    def _setEvclass(self, evclass):
        if evclass == 128:
            self._evclass = evclass
        elif evclass == 256:
            self._evclass = evclass
        elif evclass == 512:
            self._evclass = evclass
        elif evclass == 1024:
            self._evclass = evclass
        else:
            raise ValueError('Invalid Event Class. Options: Source[128] Clean[256] UltraClean[512] UltraCleanVeto[1024]')

    def getEvtype(self):
        return 3
        
    def _setEvtype(self, evtype):
        if evtype == 3:
            self._evtype = evtype
        elif evtype == 2:
            self._evtype = evtype
        elif evtype == 1:
            self._evtype = evtype
        else:
            raise ValueError('Invalid Event Type. Options: FRONT[1] BACK[2] FRONT/BACK[3] <recommended>')
       
    def getRa(self):
        return self._ra
        
    def _setRa(self, ra):
        if ra >= 0.0:
            if ra <= 360.0:
                self._ra = ra
            else:
                raise ValueError('Invalid value for right ascension of center. Valid Range: 0 to 360 degrees')
        else:
            raise ValueError('Invalid value for right ascension of center. Valid Range: 0 to 360 degrees')
            
    def getDec(self):
        return self._dec
        
    def _setDec(self,dec):
        if dec >= -90.0:
            if dec <= 90.0:
                self._dec = dec
            else:
                raise ValueError('Invalid value for declination of center. Valid Range: -90 to 90 degrees')
        else:
            raise ValueError('Invalid value for declination of center. Valid Range: -90 to 90 degrees')

    def getRad(self):
        return self._rad
        
    def _setRad(self, rad):
        if rad >= 0.0:
            if rad <= 180.0:
                self._rad = rad
            else:
                raise ValueError('Invalid value for search radius. Valid Range: 0 to 180 degrees')
        else:
            raise ValueError('Invalid value for search radius. Valid Range: 0 to 180 degrees')

    def getEmin(self):
        return self._emin
        
    def _setEmin(self, emin):
        if emin >= 30.0:
            if emin <= 2e6:
                self._emin = emin
            else:
                raise ValueError('Invalid value for minimum energy. Valid Range: 30 MeV to 2 TeV')
        else:
            raise ValueError('Invalid value for minimum energy. Valid Range: 30 MeV to 2 TeV')
        if emin <= 3e5 and emin >= 100:
            pass
        else:
            print "WARNING: Minimum energy " + str(emin/1000.0) + " GeV is outside of standard Fermi analysis range (0.1 - 300 GeV)"
        
    def getEmax(self):
        return self._emax
        
    def _setEmax(self, emax):
        if emax >= 30.0:
            if emax <= 2e6:
                self._emax = emax
            else:
                raise ValueError('Invalid value for maximum energy. Valid Range: 30 MeV to 2 TeV')
        else:
            raise ValueError('Invalid value for maximum energy. Valid Range: 30 MeV to 2 TeV')
        if emax <= 3e5 and emax >= 100:
            pass
        else:
            print "WARNING: Maximum energy " + str(emax/1000.0) + " GeV is outside of standard Fermi analysis range (0.1 - 300 GeV)"
        
    def getZmax(self):
        return self._zmax
        
    def _setZmax(self, zmax):
        if zmax >= 0.0:
            if zmax <= 180.0:
                self._zmax = zmax
                if zmax > 90:
                    print "WARNING: zmax is above recommended value of 90 degrees"
            else:
                raise ValueError('Invalid value for maximum zenith. Valid Range: 0 to 180 degrees')
        else:
            raise ValueError('Invalid value for maximum zenith. Valid Range: 0 to 180 degrees')
        
    def getTmin(self):
        return self._tmin
        
    def _setTmin(self, tmin):
        if tmin >= 239556961:
            self._tmin = tmin
        else:
            raise ValueError('Invalid value for start time. Valid Range: 239557417 to infinity MET')
        
    def getTmax(self):
        return self._tmax
        
    def _setTmax(self, tmax):
        if tmax >= 239556961:
            self._tmax = tmax
        else:
            raise ValueError('Invalid value for end time. Valid Range: 239557417 to infinity MET')
        
    def getInfile(self):
        return self._infile
        
    def _setInfile(self, infile):
        if isinstance(infile, basestring):
            self._infile = infile
        else:
            raise ValueError('Invalid python object for file name. File must be represented as string. Example: \'data.fits\'')
        
    def getOutfile(self):
        return self._outfile
        
    def _setOutfile(self, outfile):
        if isinstance(outfile, basestring):
            self._outfile = outfile
        else:
            raise ValueError('Invalid python object for file name. File must be represented as string. Example: \'data.fits\'')
    
    def getEvfile(self):
        return self._evfile
        
    def _setEvfile(self, evfile):
        if isinstance(evfile, basestring):
            self._evfile = evfile
        else:
            raise ValueError('Invalid python object for file name. File must be represented as string. Example: \'data.fits\'')
    
    def getScfile(self):
        return self._scfile
        
    def _setScfile(self, scfile):
        if isinstance(scfile, basestring):
            self._scfile = scfile
        else:
            raise ValueError('Invalid python object for file name. File must be represented as string. Example: \'spacecraft.fits\'')

    def getExpcube(self):
        return self._expcube
        
    def _setExpcube(self, expcube):
        if isinstance(expcube, basestring):
            self._expcube = expcube
        else:
            raise ValueError('Invalid python object for file name. File must be represented as string. Example: \'expcube.fits\'')

    def getDcostheta(self):
        return self._dcostheta

    def _setDcostheta(self, dcostheta):
        if dcostheta <= 1:
            if dcostheta >= 0:
                   self._dcostheta = dcostheta
        else:
            raise ValueError('Invalid vaue for python object. Costheta must be between 0 and 1.')

    def getBinsz(self):
        return self._binsz

    def _setBinsz(self, binsz):
        if binsz < 1.0:
            self._binsz = binsz
        else:
            raise ValueError('Invalid vaue for python object. Binsz must be between 0 and 1.')

    def getIrfs(self):
        return self._irfs

    def _setIrfs(self, irfs):
        if isinstance(irfs, basestring):
            self._irfs = irfs
        else:
            raise ValueError('Invalid value for python object. Irfs must be a string. Example: \'P8R2_SOURCE_V6\'')

    def getSrcrad(self):
        return self._srcrad

    def _setSrcrad(self, srcrad):
        self._srcrad = srcrad

    def getNlong(self):
        return self._nlong

    def _setNlong(self, nlong):
        if nlong <= 1000 and nlong >= 2:
            self._nlong = nlong
        else:
            raise ValueError('Invalid vaue for python object. Nlong must be between 2 and 1000.')

    def getNlat(self):
        return self._nlat

    def _setNlat(self, nlat):
        if nlat <= 1000 and nlat >= 2:
            self._nlat = nlat
        else:
            raise ValueError('Invalid vaue for python object. Nlat must be between 2 and 1000.')

    def getNenergies(self):
        return self._nenergies

    def _setNenergies(self, nenergies):
        if nenergies >= 2 and nenergies <= 100:
            self._nenergies = nenergies
        else:
            raise ValueError('Invalid value for python object. Nenergies must be between 2 and 100.')

    def getSrcmdl(self):
        return self._srcmdl

    def _setSrcmdl(self, srcmdl):
        if isinstance(srcmdl, basestring):
            self._srcmdl = srcmdl
        else:
            raise ValueError('Invalid value for python object. Srcmdl must be a sting. Example: \'yourmodel.xml\'')

    #Set to defult setting suggested by Fermi for transient
    def getFilter(self):
        return self._Filter

    def _setFilter(self, Filter):
        if isinstance(Filter, basestring):
            self._Filter = Filter
        else:
            raise ValueError('Invalid value for python object. Filter must be a sting. Example: \'(DATA_QUAL>0)&&(LAT_CONFIG==1)\'')
        
    def getRoicut(self):
        return self._roicut
        
    def _setRoicut(self, roicut):
        if isinstance(roicut, basestring):
            if roicut == 'yes':
                self._roicut = roicut
            elif roicut == 'no':
                self._roicut = roicut
        else:
            raise ValueError('Invalid entry. Must input string as \'no\' or \'yes\'.') 
   
    ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!####
    ####     All of the get and _set functions below have not been properly added to the show all      ####
    ####     function below. Nor do they conatian proper protection from user error.                   ####
    ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!####
    
    def getAlgorithm(self):
        return self._algorithm

    def _setAlgorithm(self, algorithm):
        if isinstance(algorithm, basestring):
            if algorithm == 'CCUBE':
                self._algorithm = algorithm
            elif algorithm == 'CMAP':
                self._algorithm = algorithm
            elif algorithm == 'LC':
                self._algorithm = algorithm     
            elif algorithm == 'PHA1':
                self._algorithm = algorithm     
            elif algorithm == 'PHA2':
                self._algorithm = algorithm     
            elif algorithm == 'HEALPIX':
                self._algorithm = algorithm
            else:
                raise ValueError('Invalid entry. Must input string as \'CCUBE\', \'CMAP\', \'LC\', \'PHA1\', \'PHA2\', or \'HEALPIX\'.')
        else:
            raise ValueError('Invalid entry. Must input string as \'CCUBE\', \'CMAP\', \'LC\', \'PHA1\', \'PHA2\', or \'HEALPIX\'.')

    def getCoordsys(self):
        return self._coordsys

    def _setCoordsys(self, coordsys):
        self._coordsys = coordsys

    def getNxpix(self):
        return self._nxpix

    def _setNxpix(self, nxpix):
        self._nxpix = nxpix

    def getNypix(self):
        return self._nypix

    def _setNypix(self, nypix):
        self._nypix = nypix

    def getAxisrot(self):
        return self._axisrot

    def _setAxisrot(self, axisrot):
        self._axisrot = axisrot

    def getProj(self):
        return self._proj

    def _setProj(self, proj):
        self._proj = proj

    def getEnumbins(self):
        return self._enumbins

    def _setEnumbins(self, enumbins):
        self._enumbins = enumbins

    def getCmap(self):
        return self._cmap

    def _setCmap(self, cmap):
        self._cmap = cmap

    def getBexpmap(self):
        return self._bexpmap

    def _setBexpmap(self, bexpmap):
        self._bexpmap = bexpmap
    
    def getSrcmaps(self):
        return self._srcmaps

    def _setSrcmaps(self, srcmaps):
        self._srcmaps = srcmaps

    def getOptimizer(self):
        return self._optimizer

    def _setOptimizer(self, optimizer):
        self._optimizer = optimizer

###############################################################
#          
#   FermiTools using python framework. 
#
#     Operations include:
#          gtselect
#          gtmktime
#          gtexpcube *
#          gtexpmap
#          gtdiffresps
#          gtbin(cmap/ccube)
#          gtexpcube2
#          gtsrcmaps
#          gttsmap
#          gtmodel
#   
#     *outdated - for binned analysis with new data use
#                    gtexpcube2
##############################################################
    def amonSelect(self):    
        my_apps.filter['evclass'] = self.getEvclass()
        my_apps.filter['evtype'] = self.getEvtype()
        my_apps.filter['ra'] = self.getRa()
        my_apps.filter['dec'] = self.getDec()
        my_apps.filter['rad'] = self.getRad()
        my_apps.filter['emin'] = self.getEmin()
        my_apps.filter['emax'] = self.getEmax()
        my_apps.filter['zmax'] = self.getZmax()
        my_apps.filter['tmin'] = self.getTmin()
        my_apps.filter['tmax'] = self.getTmax()
        my_apps.filter['infile'] = self.getInfile()
        my_apps.filter['outfile'] = self.getOutfile()
        my_apps.filter.run()
        
    def amonTime(self):
        my_apps.maketime['scfile'] = self.getScfile()
        my_apps.maketime['filter'] = self.getFilter()
        my_apps.maketime['roicut'] = self.getRoicut()
        my_apps.maketime['evfile'] = self.getEvfile()
        my_apps.maketime['outfile'] = self.getOutfile()
        my_apps.maketime.run()

    def amonExpcube(self):
        my_apps.expCube['evfile'] = self.getEvfile()
        my_apps.expCube['scfile'] = self.getScfile()
        my_apps.expCube['outfile'] = self.getOutfile()
        my_apps.expCube['zmax'] = self.getZmax()
        my_apps.expCube['dcostheta'] = self.getDcostheta()
        my_apps.expCube['binsz'] = self.getBinsz()
        my_apps.expCube.run()

    def amonExpmap(self):
        my_apps.expMap['evfile'] = self.getEvfile()
        my_apps.expMap['scfile'] = self.getScfile()
        my_apps.expMap['expcube'] = self.getExpcube()
        my_apps.expMap['outfile'] = self.getOutfile()
        my_apps.expMap['irfs'] = self.getIrfs()
        my_apps.expMap['srcrad'] = self.getSrcrad()
        my_apps.expMap['nlong'] = self.getNlong()
        my_apps.expMap['nlat'] = self.getNlat()
        my_apps.expMap['nenergies'] = self.Nenergies()
        my_apps.expMap.run()

    def amonDiffresps(self):
        my_apps.diffResps['evfile'] = self.getEvfile()
        my_apps.diffResps['scfile'] = self.getScfile()
        my_apps.diffResps['srcmdl'] = self.getSrcmdl()
        my_apps.diffResps['irfs'] = self.getIrfs()
        my_apps.diffResps.run()

    def amonBincmap(self):
        my_apps.evtbin['algorithm'] = self.getAlgorithm()
        my_apps.evtbin['evfile'] = self.getEvfile()
        my_apps.evtbin['outfile'] = self.getOutfile()
        my_apps.evtbin['scfile'] = self.getScfile()
        my_apps.evtbin['nxpix'] = self.getNxpix()
        my_apps.evtbin['nypix'] = self.getNypix()
        my_apps.evtbin['binsz'] = self.getBinsz()
        my_apps.evtbin['coordsys'] = self.getCoordsys()
        my_apps.evtbin['xref'] = self.getRa()
        my_apps.evtbin['yref'] = self.getDec()
        my_apps.evtbin['axisrot'] = self.getAxisrot()
        my_apps.evtbin['proj'] = self.getProj()
        my_apps.evtbin['emin'] = self.getEmin()
        my_apps.evtbin['emax'] = self.getEmax()
        my_apps.evtbin.run()

    def amonBinccube(self):
        my_apps.evtbin['algorithm'] = self.getAlgorithm()
        my_apps.evtbin['evfile'] = self.getEvfile()
        my_apps.evtbin['outfile'] = self.getOutfile()
        my_apps.evtbin['scfile'] = self.getScfile()
        my_apps.evtbin['nxpix'] = self.getNxpix()
        my_apps.evtbin['nypix'] = self.getNypix()
        my_apps.evtbin['binsz'] = self.getBinsz()
        my_apps.evtbin['coordsys'] = self.getCoordsys()
        my_apps.evtbin['xref'] = self.getRa()
        my_apps.evtbin['yref'] = self.getDec()
        my_apps.evtbin['axisrot'] = self.getAxisrot()
        my_apps.evtbin['proj'] = self.getProj()
        my_apps.evtbin['ebinalg'] = 'LOG'
        my_apps.evtbin['enumbins'] = self.getEnumbins()
        my_apps.evtbin['emin'] = self.getEmin()
        my_apps.evtbin['emax'] = self.getEmax()
        my_apps.evtbin.run()

    def amonExpcube2(self):
        my_apps.gtexpcube2['infile'] = self.getInfile()
        my_apps.gtexpcube2['cmap'] = 'none'
        my_apps.gtexpcube2['nxpix']= self.getNxpix()
        my_apps.gtexpcube2['nypix'] = self.getNypix()
        my_apps.gtexpcube2['binsz'] = self.getBinsz()
        my_apps.gtexpcube2['coordsys'] = self.getCoordsys()
        my_apps.gtexpcube2['xref'] = self.getRa()
        my_apps.gtexpcube2['yref'] = self.getDec()
        my_apps.gtexpcube2['axisrot'] = self.getAxisrot()
        my_apps.gtexpcube2['proj'] = self.getProj()
        my_apps.gtexpcube2['emin'] = self.getEmin()
        my_apps.gtexpcube2['emax'] = self.getEmax()
        my_apps.gtexpcube2['enumbins'] = self.getEnumbins()
        my_apps.gtexpcube2['outfile'] = self.getOutfile()
        my_apps.gtexpcube2['irfs'] = self.getIrfs()
        my_apps.gtexpcube2['clobber'] = 'yes'
        my_apps.gtexpcube2.run()

    def amonSrcmaps(self):
        my_apps.srcMaps['scfile'] = self.getScfile()
        my_apps.srcMaps['expcube'] = self.getExpcube()
        my_apps.srcMaps['cmap'] = self.getCmap()
        my_apps.srcMaps['bexpmap'] = self.getBexpmap()
        my_apps.srcMaps['srcmdl'] = self.getSrcmdl()
        my_apps.srcMaps['outfile'] = self.getOutfile()
        my_apps.srcMaps['irfs'] = self.getIrfs()
        my_apps.srcMaps.run()

    def amonModelmaps(self):
        my_apps.model_map['srcmaps'] = self.getSrcmaps()
        my_apps.model_map['expcube'] = self.getExpcube()
        my_apps.model_map['bexpmap'] = self.getBexpmap()
        my_apps.model_map['srcmdl'] = self.getSrcmdl()
        my_apps.model_map['outfile'] = self.getOutfile()
        my_apps.model_map['irfs'] = self.getIrfs()
        my_apps.model_map['evtype'] = self.getEvtype()
        my_apps.model_map.run()

    def amonTsmap(self):
        my_apps.TsMap['statistic'] = "BINNED"
        my_apps.TsMap['scfile'] = self.getScfile()
        my_apps.TsMap['cmap'] = self.getCmap()
        my_apps.TsMap['evfile'] = self.getEvfile()
        my_apps.TsMap['bexpmap'] = self.getBexpmap()
        my_apps.TsMap['expcube'] = self.getExpcube()
        my_apps.TsMap['srcmdl'] = self.getSrcmdl()
        my_apps.TsMap['irfs'] = self.getIrfs()
        my_apps.TsMap['optimizer'] = self.getOptimizer()
        my_apps.TsMap['outfile'] = self.getOutfile()
        my_apps.TsMap['nxpix'] = self.getNxpix()
        my_apps.TsMap['nypix'] = self.getNypix()
        my_apps.TsMap['binsz'] = self.getBinsz()
        my_apps.TsMap['coordsys'] = self.getCoordsys()
        my_apps.TsMap['xref'] = self.getRa()
        my_apps.TsMap['yref'] = self.getDec()
        my_apps.TsMap['proj'] = 'AIT'
        my_apps.TsMap.run()

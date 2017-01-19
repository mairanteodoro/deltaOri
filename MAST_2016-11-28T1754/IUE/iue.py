#!/usr/bin/env python

class Iue():

    # list with the filename (full path) for all data to be used.
    # (class attribute; will be changed to tuple by .readIueData())
    myList = ''

    def __init__(self, *arg):
        # arg == myList
        # N.B.: the Iue class can be instantiated with or without providing a list of the files to be read and stored
        try:
            if (len(arg) == 1):
                # if arg has been provided then read data immediately
                self.readIueData(arg[0])
            else:
                return
        except:
            pass

    def __checkMyList__(self):
        '''
            This method checks for the existence of IUE data.
        '''
        import sys

        try:
            if (type(self.myList) != tuple) or (len(self.myList[0]) == 0):
                print("")
                print("### ")
                print("### WARNING: No IUE data has been loaded yet.")
                print("###")
                print("### Please, make sure to run {0}.readIueData(list) first!".format(self.__class__.__name__))
                print("### ")
                print("")
                raise SystemExit()
        except SystemExit:
            sys.exit()

    def readIueData(self, myList):
        '''
            This method appends 2 attributes to self: self.myInfo and self.mySpec. They correspond to the info about the observation/data and the full spectrum, respectively. Both are saved in pandas data frame format.

            - In order to access the information about the files:
                self.myInfo.head()

            - In order to access the full spectrum of a given observation:
                - first find the index of the desired file:
                idx = self.myInfo[self.myInfo.filename.str.contains("swp42749")].index[0]
                - now retrieve the wavelength and flux arrays:
                self.mySpec[idx].wave and self.mySpec[idx].flux
        '''

        import pandas as pd
        import astropy.units as u
        from astropy.time import Time, TimeDelta
        import astropy.io.fits as fits
        import numpy as np

        self.myList = myList, # <- from now on, this is a tuple!

        dataToBeStored = [
            # [filename, PrimaryHDU, BinTableHDU]
            myList,
            # [ list(fits.open(x))[0] for x in myList ],
            # [ list(fits.open(x))[1] for x in myList ],
            # info about observation
            #
            # APERTURE
            [ fits.open(x)[0].header.get('APERTURE') for x in myList ],
            # DATE + TIME
            [ fits.open(x)[0].header.get('SDATEOBS' if fits.open(x)[0].header.get('SDATEOBS') else 'LDATEOBS') + 'T' + fits.open(x)[0].header.get('STIMEOBS' if fits.open(x)[0].header.get('STIMEOBS') else 'LTIMEOBS') for x in myList ],
            # EXPTIME
            [ fits.open(x)[0].header.get('SEXPTIME' if fits.open(x)[0].header.get('SEXPTIME') else 'LEXPTIME') for x in myList ],
            # JD
            [ fits.open(x)[0].header.get('SJD-OBS' if fits.open(x)[0].header.get('SJD-OBS') else 'LJD-OBS') for x in myList ],
            # MID JD
            [ fits.open(x)[0].header.get('SJD-MID' if fits.open(x)[0].header.get('SJD-MID') else 'LJD-MID') for x in myList ]
        ]

        self.myInfo = pd.DataFrame(data=dataToBeStored)
        # transpose
        self.myInfo = self.myInfo.transpose()
        # set column names
        self.myInfo.columns = [
            "filename",
            # "hdu",
            # "binTable",
            # info about observation
            "aperture",
            "DateObs",
            "exptime",
            "jdObs",
            "jdMid",
            ]
        # converting DATE + TIME into datetime type so that it can be handled like self.myInfo.datetime.year
        self.myInfo['DateObs'] = pd.to_datetime(self.myInfo['DateObs'], infer_datetime_format=True)
        # read each FITS file and store the corresponding full spectrum
        self.mySpec = [ self.getSpec(x) for x in myList ]

    def getSpec(self, filename):
        '''
            This method loops over the orders of a given file and returns a dataframe with columns 'order', 'wave', 'flux', and 'filename'.
        '''

        self.__checkMyList__()

        from astropy.table import Table
        import pandas as pd
        import numpy as np
        import astropy.units as u

        # find the filename's index
        fIndex = self.findIndexOf(filename)
        # read and return data in tabular format
        tb = Table.read(self.myInfo.filename[fIndex], hdu=1)
        wave, flux = [], []

        order = [tb['ORDER'][i] for i in range(tb['ORDER'].shape[0])]
        for index in range(len(order)):
            # wavelength
            wave.append(tb['WAVELENGTH'][index] + ( range((tb['RIPPLE'][index]).shape[0]) - tb['STARTPIX'][index] ) * tb['DELTAW'][index])
            # flux
            flux.append(tb['RIPPLE'][index])

        '''
        Return a dataframe with the spectrum for each order:
             order   wave    flux    filename
             -----   -----   -----   -----
        0    #       []      []      string
        1    #       []      []      string
        '''
        return pd.DataFrame({'filename': filename,
                            'order': order,
                            'wave': wave,
                            'flux': flux})

    def getAllSpec (self):
        '''
            Return an array where each element is a pandas DataFrame containing data about 'order', 'wave', and 'flux' for the corresponding 'filename'.

            e.g.:
                iue = Iue()
                iue.readIueData(file)
                allSpec = iue.getAllSpec()
                allSpec[0].head()
                allSpec[0].filename
                allSpec[0].order, allSpec[0].wave, allSpec[0].flux
        '''

        self.__checkMyList__()

        # return an array where each element is a pandas DataFrame
        # containing 'filename', 'order', 'wave', and 'flux'
        return [ self.getSpec(filename) for filename in self.myList[0] ]

    def __mergeOrders__(self):
        '''
            Loop over the orders and splice them together where they overlap.
        '''

        self.__checkMyList__()

        # loop over the files
        for item in self.myList[0]:
            # get the spectral data corresponding to the current file
            currSpec = self.mySpec[self.findIndexOf(item)]
            # loop over the orders
            for order in currSpec.order:
                # do the splicing here
                print(order)

    def findIndexOf(self, pattern):
        '''
            Find the unique index corresponding to the name of the file that matches up the provided pattern.

            e.g.: index = iue.findIndexOf("05722")
        '''

        self.__checkMyList__()

        import numpy as np

        try:
            return self.myInfo[self.myInfo.filename.str.contains(pattern)].index[0]
        except:
            print("findIndexOf exception: no match found.")
            # return None (type(returnedValue) = None)
            pass

    def findByDate(self, start, end):
        '''
            Search dataframe for data within a given interval in date format.
        '''

        self.__checkMyList__()

        s = start if start < end else end
        e = end if start < end else start
        mask = (self.myInfo['DateObs'] >= s) & (self.myInfo['DateObs'] <= e)
        return self.myInfo.loc[mask]


class IueBinSys(Iue):
    '''
        This class will handle data from binary systems by adding an extra column containing the orbital phase of the observation.

        Usage:
            binSysName = IueBinSys() # new instance of the class
            binSysName.readIueData(fileList) # read data list
            binSysName.setEphemeris(period, T0) # set ephemeris for the system
            binSysName.myInfo.head() # show the resulting dataframe
    '''

    def __init__(self):
        pass

    def setEphemeris(self, period, T0):
        '''
            Set the period and zero point for orbital phase calculation.
        '''

        self.checkMyList()

        import astropy.units as u
        from astropy.time import Time, TimeDelta

        self.per = TimeDelta(period * u.day)
        self.T0 = Time(T0, format='jd')
        # ephemeris
        # T - T0
        deltaT = self.myInfo['jdMid'] - self.T0.jd
        self.myInfo['deltaT'] = deltaT
        # phase = mod( (T-T0) / per )           V  (T-T0) > 0
        # phase = mod(mod( (T-T0) / per ) + 1)  V  (T-T0) < 0
        self.myInfo['phase'] = [ ( (dT/self.per.value)%1 if (dT > self.T0.jd) else ((dT/self.per.value)%1+1)%1 ) for dT in deltaT]

    def findByPhase(self, start, end):
        '''
            Search dataset for data within a given interval in phase.
        '''

        self.checkMyList()

        s = start if start < end else end
        e = end if start < end else start
        mask = (self.myInfo['phase'] >= s) & (self.myInfo['phase'] <= e)
        return self.myInfo.loc[mask]


if __name__ == "__main__":

    print("### Usage ###")
    print("")
    print("For general IUE data:")
    print(" -> form iue import Iue")
    print(" -> obj = Iue()")
    print("")
    print("For observation of binary systems:")
    print(" -> from iue import IueBinSys")
    print(" -> obj = IueBinSys()")

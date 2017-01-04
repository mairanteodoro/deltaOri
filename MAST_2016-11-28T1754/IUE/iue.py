#!/usr/bin/env python

class Iue():

    def __init__(self, myList):
        import pandas as pd
        import astropy.units as u
        from astropy.time import Time, TimeDelta
        import astropy.io.fits as fits
        import numpy as np

        dataToBeStored = [
            # [filename, PrimaryHDU, BinTableHDU]
            myList,
            [ list(fits.open(x))[0] for x in myList ],
            [ list(fits.open(x))[1] for x in myList ],
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
            "hdu",
            "binTable",
            # info about observation
            "aperture",
            "myDatetime",
            "exptime",
            "jdObs",
            "jdMid",
            ]
        # converting DATE + TIME into datetime type so that it can be handled like self.myInfo.datetime.year
        self.myInfo['myDatetime'] = pd.to_datetime(self.myInfo['myDatetime'], infer_datetime_format=True)
        # ephemeris from Corcoran et al. 2015
        per = TimeDelta(5.732436 * u.day)
        T0 = Time(2456295.674, format='jd')
        # T - T0
        deltaT = self.myInfo['jdMid'] - T0.jd
        self.myInfo['deltaT'] = deltaT
        # phase = mod( (T-T0) / per )           V  (T-T0) > 0
        # phase = mod(mod( (T-T0) / per ) + 1)  V  (T-T0) < 0
        self.myInfo['phase'] = [ ( (dT/per.value)%1 if (dT > T0.jd) else ((dT/per.value)%1+1)%1 ) for dT in deltaT]
        # read each FITS file and store all the full spectrum
        self.mySpec = [ self.createSpec(x) for x in myList ]

    def createSpec(self, filename):
        '''
            This method loops over the orders
            of a given file and return the full spectrum.
        '''
        from astropy.table import Table
        import pandas as pd
        import numpy as np
        import astropy.units as u

        # find the filename's index
        fIndex = self.myInfo.filename[
            self.myInfo.filename == filename
        ].index[0]
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
        Return a dataframe with the
        spectrum for each order:
             order   wave    flux
             -----   -----   -----
        0    #       []      []
        1    #       []      []
        '''
        return pd.DataFrame({'filename': filename,
                            'order': order,
                            'wave': wave,
                            'flux': flux})

    def getAllSpec (self, fileList):
        '''
            Return an array where each element is a pandas DataFrame containing data about 'order', 'wave', and 'flux' for the corresponding 'filename'.

            e.g.:
                allSpec = iue.getAllSpec("./iueData.json")
                allSpec[0].head()
                allSpec[0].filename
                allSpec[0].order, allSpec[0].wave, allSpec[0].flux
        '''
        import json as json

        # test if fileList is of type dict or list;
        # if not, load it into a dict and select the 'filename' key
        myList = fileList if (type(fileList) == dict or type(fileList) == list) else (json.load( open(fileList) ))['filename']

        # return an array where each element is a pandas DataFrame
        # containing 'filename', 'order', 'wave', and 'flux'
        return [ self.createSpec(filename) for filename in myList ]



###
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
import json as json

# parse filenames into list
fileList = json.load( open( "./iueData.json" ) )
# instantiate an Iue object containing
# all data inside the instance attribute 'myData'
iue = Iue( fileList['filename'] )

# # read and store all spectra
# # using a JSON file
# allSpec = iue.getAllSpec("./iueData.json")
# # using the output from json.load()
# allSpec2 = iue.getAllSpec(fileList['filename'])




# get the distribution of phase coverage
# plt.ion()
# iue.myData.hist(column='phase')




# get all the orders for a file
# spec = iue.createSpec("./DATA/lwr02241.mxhi")

# # plot all the orders
# fig, ax = plt.subplots()
# [ax.plot(res['wave'][i], res['flux'][i], marker='') for i in range(len(res['wave']))]
# fig.savefig('./RESULTS/spec.pdf')
# plt.clf()
# plt.close('all')
#
# # plt.ion()
# # # loop over the orders and plot them
# # for index in range(61):
# #     plt.plot(res[0][index], res[1][index], marker='')

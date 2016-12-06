#!/usr/bin/env python

class Iue():

    def __init__(self, myList):
        import pandas as pd
        import astropy.io.fits as fits

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
            [ fits.open(x)[0].header.get('SJD-MID') if fits.open(x)[0].header.get('SJD-MID') else 'LJD-MID' for x in myList ],
        ]

        self.myData = pd.DataFrame(data=dataToBeStored)
        # transpose
        self.myData = self.myData.transpose()
        # set column names
        self.myData.columns = [
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
        # converting DATE + TIME into datetime
        # type that can be handled like self.myData.datetime.year
        self.myData['myDatetime'] = pd.to_datetime(self.myData['myDatetime'], infer_datetime_format=True)

    def createSpec(self, filename):
        '''
            Ideally, this method should create
            a full spectrum for a given filename.
            This method should loop over the orders
            of a given file and return the full spectrum.
        '''
        from astropy.table import Table
        import pandas as pd
        import numpy as np
        # find the filename's index
        fIndex = self.myData.filename[
            self.myData.filename == filename
        ].index[0]
        # read and return data in tabular format
        tb = Table.read(self.myData.filename[fIndex], hdu=1)
        '''
            Instead of a list:
            - create a pandas dataframe?
            - containing what?
        '''
        wave, flux = [], []
        for index in range(61):
            # wavelength
            wave.append(tb['WAVELENGTH'][index] + ( range((tb['RIPPLE'][index]).shape[0]) - tb['STARTPIX'][index] ) * tb['DELTAW'][index])
            # flux
            flux.append(tb['RIPPLE'][index])

        # return a dict with wavelength and flux
        # for each echelle order (final dimension: 61, 768)
        # (e.g. len(res['wave']) = 61; len(res['wave'][0]) = 768)
        return {'wave': wave, 'flux': flux}


###
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
import json as json

# parse filenames into list
fileList = json.load( open( "./iueData.json" ) )
# instantiate an Iue object containing all data
iue = Iue( fileList['filename'] )

# get all the orders
res = iue.createSpec("./DATA/lwr02241.mxhi")
# plot all the orders
fig, ax = plt.subplots()
[ax.plot(res['wave'][i], res['flux'][i], marker='') for i in range(len(res['wave']))]
fig.savefig('./RESULTS/spec.pdf')
plt.clf()
plt.close('all')

# plt.ion()
# # loop over the orders and plot them
# for index in range(61):
#     plt.plot(res[0][index], res[1][index], marker='')

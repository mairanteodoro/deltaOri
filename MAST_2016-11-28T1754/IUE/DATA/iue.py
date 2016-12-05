#!/usr/bin/env python

class Iue():

    def __init__(self, myList):
        import pandas as pd
        import astropy.io.fits as fits

        # [filename, PrimaryHDU, BinTableHDU]
        self.data = pd.DataFrame([
            myList,
            [ list(fits.open(x))[0] for x in myList ],
            [ list(fits.open(x))[1] for x in myList ]
        ])
        # transpose
        self.data = self.data.transpose()
        # set column names
        self.data.columns = ["filename", "hdu", "binTable"]

    def createSpec(self, filename):
        '''
            Ideally, this method should create
            a full spectrum for a given filename.
            This method should loop over the orders
            of a given file and return the full spectrum.
            (to merge or not to merge?)
        '''
        from astropy.table import Table
        # find the filename's index
        fIndex = self.data.filename[
            self.data.filename == filename
        ].index[0]
        # read and return data in tabular format
        tb = Table.read(self.data.filename[fIndex], hdu=1)
        '''
            Instead of a list, this method should
            create a pandas dataframe containing what?
        '''
        wave, flux = [], []
        for index in range(61):
            # wavelength
            wave.append(tb['WAVELENGTH'][index] + ( range((tb['RIPPLE'][index]).shape[0]) - tb['STARTPIX'][index] ) * tb['DELTAW'][index])
            # flux
            flux.append(tb['RIPPLE'][index])
        # return list with wavelength and flux
        return [wave, flux]



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

# get spectrum from file
res = iue.createSpec("./DATA/lwr02241.mxhi")

plt.ion()
# loop over the orders and plot them
for index in range(61):
    plt.plot(res[0][index], res[1][index], marker='')

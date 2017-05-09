#!/usr/bin/env python

class Iue(object): # <- providing 'object' to comply with new style class!

    # flag for the existence of list; class attribute
    _my_list_ok = False

    def __init__(self, *arg):
        # arg == my_list
        # N.B.: the Iue class can be instantiated with or without providing a list of the files to be read in and stored
        try:
            if (arg.__len__() == 1):
                # if arg has been provided then read data immediately
                self.read_iue_data(arg[0])
            else:
                return
        except:
            pass

    def _check_my_list(self):
        '''
            This method checks for the existence of IUE data.
        '''
        import sys

        try:
            if (type(self.my_list) != tuple) or (len(self.my_list[0]) == 0):
                print("")
                print("### ")
                print("### WARNING: No IUE data has been loaded yet.")
                print("###")
                print("### Please, make sure to run {0}.read_iue_data(list) first!".format(self.__class__.__name__))
                print("### ")
                print("")
                raise SystemExit()
        except SystemExit:
            sys.exit()

    def read_iue_data(self, my_list):
        '''
            This method appends 2 attributes to self: self.my_info and self.my_spec. They correspond to the info about the observation/data and the full spectrum, respectively. Both are saved in pandas data frame format.

            - In order to access the information about the files:
                self.my_info.head()

            - In order to access the full spectrum of a given observation:
                - first find the index of the desired file:
                idx = self.my_info[self.my_info.filename.str.contains("swp42749")].index[0]
                - now retrieve the wavelength and flux arrays:
                self.my_spec[idx].wave and self.my_spec[idx].flux
        '''

        import pandas as pd
        import astropy.units as u
        from astropy.time import Time, TimeDelta
        import astropy.io.fits as fits
        import numpy as np

        self.my_list = my_list, # <- this is a tuple!

        data_to_be_stored = [
            # [filename, PrimaryHDU, BinTableHDU]
            my_list,
            # [ list(fits.open(x))[0] for x in my_list ],
            # [ list(fits.open(x))[1] for x in my_list ],
            # info about observation
            #
            # APERTURE
            [ fits.open(x)[0].header.get('APERTURE') for x in my_list ],
            # DATE + TIME
            [ fits.open(x)[0].header.get('SDATEOBS' if fits.open(x)[0].header.get('SDATEOBS') else 'LDATEOBS') + 'T' + fits.open(x)[0].header.get('STIMEOBS' if fits.open(x)[0].header.get('STIMEOBS') else 'LTIMEOBS') for x in my_list ],
            # EXPTIME
            [ fits.open(x)[0].header.get('SEXPTIME' if fits.open(x)[0].header.get('SEXPTIME') else 'LEXPTIME') for x in my_list ],
            # JD
            [ fits.open(x)[0].header.get('SJD-OBS' if fits.open(x)[0].header.get('SJD-OBS') else 'LJD-OBS') for x in my_list ],
            # MID JD
            [ fits.open(x)[0].header.get('SJD-MID' if fits.open(x)[0].header.get('SJD-MID') else 'LJD-MID') for x in my_list ]
        ]

        self.my_info = pd.DataFrame(data=data_to_be_stored)
        # transpose
        self.my_info = self.my_info.transpose()
        # set column names
        self.my_info.columns = [
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
        # converting DATE + TIME into datetime type so that it can be handled like self.my_info.datetime.year
        self.my_info['DateObs'] = pd.to_datetime(self.my_info['DateObs'], infer_datetime_format=True)
        # set flag to True
        self._my_list_ok = True
        # read each FITS file and store the corresponding full spectrum
        self.my_spec = [ self.get_spec(x) for x in my_list ]

    def get_spec(self, filename):
        '''
            This method loops over the orders of a given file.

            Returns a dataframe with columns 'order', 'wave', 'flux', and 'filename'.
        '''

        if (self._my_list_ok == False):
            self._check_my_list()

        from astropy.table import Table
        import pandas as pd
        import numpy as np
        import astropy.units as u

        # find the filename's index
        fIndex = self.find_index_of(filename)
        # read and return data in tabular format
        tb = Table.read(self.my_info.filename[fIndex], hdu=1)
        wave, ripple, abs_cal, wave_min, wave_max, nu = [], [], [], [], [], []

        order = [tb['ORDER'][i] for i in range(tb['ORDER'].shape[0])]
        for index in range(len(order)):
            # wavelength
            wave.append(tb['WAVELENGTH'][index] + ( range((tb['RIPPLE'][index]).shape[0]) - tb['STARTPIX'][index] ) * tb['DELTAW'][index])
            # wavelength intervals
            wave_min.append(wave[-1][0])
            wave_max.append(wave[-1][-1])
            # flux
            ripple.append(tb['RIPPLE'][index])
            abs_cal.append(tb['ABS_CAL'][index])
            nu.append(tb['QUALITY'][index])

        '''
        Return a dataframe with the spectrum for each order:
             order   wave    fluxes    filename
             -----   -----   -----   -----
        0    #       []      []      string
        1    #       []      []      string
        '''
        return pd.DataFrame({'filename': filename,
                            'order': order,
                            'wave': wave,
                            'ripple': ripple,
                            'abs_cal': abs_cal,
                            'waveMin': wave_min,
                            'waveMax': wave_max,
                            'quality': nu})

    def get_all_spec(self):
        '''
            Return an array where each element is a pandas DataFrame containing data about 'order', 'wave', and 'flux' for the corresponding 'filename'.

            e.g.:
                iue = Iue()
                iue.read_iue_data(file)
                allSpec = iue.get_all_spec()
                allSpec[0].head()
                allSpec[0].filename
                allSpec[0].order, allSpec[0].wave, allSpec[0].flux
        '''

        if (self._my_list_ok == False):
            self._check_my_list()

        # return an array where each element is a pandas DataFrame
        # containing 'filename', 'order', 'wave', and 'flux'
        return [ self.get_spec(filename) for filename in self.my_list[0] ]

    @staticmethod
    def _splice(w1, w2, f1, f2, q1, q2):
        '''
            Splice two spectra together taking into account data quality flag.

            Adapted from:
            https://archive.stsci.edu/pub/iue/software/iuedac/procedures/mergeav.pro

            *** N.B.: w1[0] MUST BE < w2[0] ***

            E.g.: w1 = iue.my_spec[0].wave.iloc[0]
                  w2 = iue.my_spec[0].wave.iloc[1]
                  f1 = iue.my_spec[0].abs_cal.iloc[0]
                  f2 = iue.my_spec[0].abs_cal.iloc[1]
                  q1 = iue.my_spec[0].quality.iloc[0]
                  q2 = iue.my_spec[0].quality.iloc[1]
                  merged_spec = iue._splice(w1, w2, f1, f2, q1, q2)
        '''
        from scipy import interpolate
        import numpy as np

        # considering only good quality flags
        good1, good2 = np.where(q1 >= 0)[0], np.where(q2 >= 0)[0]
        w1, f1, q1 = w1[good1], f1[good1], q1[good1]
        w2, f2, q2 = w2[good2], f2[good2], q2[good2]

        # OVERLAPPING REGION
        # finding the overlapping region
        w1_overlap, f1_overlap, q1_overlap = w1[w1 >= w2[0]], f1[w1 >= w2[0]], q1[w1 >= w2[0]]
        w2_overlap, f2_overlap, q2_overlap = w2[w2 <= w1[-1]], f2[w2 <= w1[-1]], q2[w2 <= w1[-1]]

        # NON OVERLAPPING REGION
        # indexes of the non-overlapping region
        non_overlap_indx = np.where(w1 >= w2[0])[0][0]
        # wave at the non-overlapping region
        w1_non_overlap, f1_non_overlap, q1_non_overlap = w1[:non_overlap_indx], f1[:non_overlap_indx], q1[:non_overlap_indx]
        w2_non_overlap, f2_non_overlap, q2_non_overlap = w2[non_overlap_indx:], f2[non_overlap_indx:], q2[non_overlap_indx:]

        # INTERPOLATION
        # determine the representation of (w2_overlap, f2_overlap)
        # using a cubic spline without smoothing (s=0)
        tck = interpolate.splrep(w2_overlap, f2_overlap, s=0)
        # interpolate the values of f2_overlap using the knots of w1_overlap
        f2new_overlap = interpolate.splev(w1_overlap, tck, der=0)

        # SPLICING
        # take the average of f1_overlap and f2new_overlap
        f_overlap = np.average([f1_overlap, f2new_overlap], axis=0)

        # OUTPUT
        # find the wavelengths of w2 that can be used to fill in the wavelengths of w1 that have q1 < 0
        w2_good2 = w2[np.where(w2>=w1[-1])[0]]
        f2_good2 = f2[np.where(w2>=w1[-1])[0]]
        # use w1_overlap because f2_overlap was interpolated onto the w1_overlap array
        w3 = np.concatenate([w1_non_overlap, w1_overlap, w2_good2])
        f3 = np.concatenate([f1_non_overlap, f_overlap, f2_good2])

        return [w3, f3]

    def _merge_orders(self):
        '''
            Loop over the orders and splice them together where they overlap.
        '''
        import time

        if (self._my_list_ok == False):
            self._check_my_list()

        # loop over the files
        start = time.time()
        for item in self.my_list[0]:
            # get the spectral data corresponding to the current file
            curr_spec = self.my_spec[self.find_index_of(item)]
            # loop over the orders
            i = 0
            for order in curr_spec.order:
                # do the splicing here
                w1, w2 = curr_spec.wave[0+i:2+i:1].values
                f1, f2 = curr_spec.flux[0+i:2+i:1].values
                res = self._splice(w1, w2, f1, f2)
                # print(curr_spec.order[0+i], curr_spec.order[2+i])
                # tempSpec1 = curr_spec.flux[0+i:2+i:1]
                i+=1
        print(time.time() - start)

    def find_index_of(self, pattern):
        '''
            Search .my_info for the index corresponding to the name of the file that matches up the provided pattern. The returned index can be used with anything within the .my_info dataframe.

            e.g.: index = iue.find_index_of("05722")
                  iue.exptime[index] # -> find exptime
                  iue.get_spec(iue.my_info.filename[index]) # -> get spectrum
        '''

        if (self._my_list_ok == False):
            self._check_my_list()

        import numpy as np

        try:
            return self.my_info[self.my_info.filename.str.contains(pattern)].index[0]
        except:
            print("find_index_of exception: no match found.")
            # return None (type(returnedValue) = None)
            pass

    def find_by_date(self, start, end):
        '''
            Search dataframe for data within a given interval in date format.

            Returns a list with a pandas dataframe for each element.
        '''
        import pandas as pd
        import numpy as np

        if (self._my_list_ok == False):
            self._check_my_list()

        s = start if start < end else end
        e = end if start < end else start

        # create a True/False mask
        mask = [(self.my_info.DateObs[i] >= pd.Timestamp(s)) & (self.my_info.DateObs[i] <= pd.Timestamp(e)) for i in range(len(self.my_info.DateObs))]

        # create list with indexes where mask == True
        mask_index = np.where(np.asarray(mask) == True)[0].tolist()

        # create series to be concatenated
        norders = [self.my_spec[mask_index[i]].shape[0] for i in range(len(mask_index))] # NORDERS IS VARIABLE!!
        nelements = len([item for item in self.my_info.jdObs.iloc[mask].values])
        # create a list with a pandas Series for each element
        s1 = [pd.Series(
            [self.my_info.DateObs.iloc[mask].values[i]] * norders[i],
            name='DateObs') for i in range(nelements)]
        s2 = [pd.Series(
            [self.my_info.aperture.iloc[mask].values[i]] * norders[i],
            name='aperture') for i in range(nelements)]
        s3 = [pd.Series(
            [self.my_info.jdObs.iloc[mask].values[i]] * norders[i],
            name='jdObs') for i in range(nelements)]
        s4 = [pd.Series(
            [self.my_info.jdMid.iloc[mask].values[i]] * norders[i],
            name='jdMid') for i in range(nelements)]
        s5 = [pd.Series(
            [self.my_info.exptime.iloc[mask].values[i]] * norders[i],
            name='exptime') for i in range(nelements)]

        # create a list with one valid data frame per element
        temp = [pd.concat([self.my_spec[mask_index[i]], s1[i], s2[i], s3[i], s4[i], s5[i]], axis=1) for i in range(len(mask_index))]

        # test for the existence of the key 'phase'
        if ('phase' in self.my_info):
            phase = [pd.Series(
                [self.my_info.phase.iloc[mask].values[i]] * norders[i], name='phase') for i in range(nelements)]
            # if present, concatenate it to the dataframe
            temp2 = [pd.concat([temp[i], phase[i]], axis=1) for i in range(len(temp))]
        else:
            # ignore it
            temp2 = temp

        # remove empty data
        res = []
        for item in temp2:
            if item.shape[0] != 0:
                res.append(item)

        # return a list with data that encompass the provided wavelength
        return res

    def find_by_wave(self, wave):
        '''
            Find data that contains the provided wavelength.

            Returns a list with a pandas dataframe for each element.
        '''
        import pandas as pd

        if (self._my_list_ok == False):
            self._check_my_list()

        # create a True/False mask
        mask = [(self.my_spec[i].waveMin < wave) & (self.my_spec[i].waveMax > wave) for i in range(len(self.my_spec))]

        # create a list with one valid data frame per element
        temp = [pd.concat( [self.my_spec[i][mask[i]], pd.Series( [self.my_info.DateObs[i] for j in range(len(self.my_spec[i][mask[i]]))], name='DateObs', index=self.my_spec[i][mask[i]].index ), pd.Series( [self.my_info.aperture[i] for j in range(len(self.my_spec[i][mask[i]]))], name='aperture', index=self.my_spec[i][mask[i]].index ), pd.Series( [self.my_info.jdObs[i] for j in range(len(self.my_spec[i][mask[i]]))], name='jdObs', index=self.my_spec[i][mask[i]].index ), pd.Series( [self.my_info.jdMid[i] for j in range(len(self.my_spec[i][mask[i]]))], name='jdMid', index=self.my_spec[i][mask[i]].index ), pd.Series( [self.my_info.exptime[i] for j in range(len(self.my_spec[i][mask[i]]))], name='exptime', index=self.my_spec[i][mask[i]].index )], axis=1 ) for i in range(len(self.my_spec))]

        # test for the existence of the key 'phase'
        if ('phase' in self.my_info):
            # if present, concatenate it to the dataframe
            temp2 = [pd.concat( [temp[i], pd.Series( [self.my_info.phase[i] for j in range(len(self.my_spec[i][mask[i]]))], name='phase', index=self.my_spec[i][mask[i]].index )], axis=1 ) for i in range(len(self.my_spec))]
        else:
            # ignore it
            temp2 = temp

        # remove empty data
        res = []
        for item in temp2:
            if item.shape[0] != 0:
                res.append(item)

        # return a list with data that encompass the provided wavelength
        return res

    def smTVS(self, *args):
        '''
        This method uses the smoothed temporal variance (smTVS) spectrum to analize light phase variations (LPV).
        
        Reference:
            (1) "Smoothed Temporal Variance Spectrum: weak line profile variations and NRP diagnostics" by Kholtygin, A. F., Sudnik, N. P. (2016)
            (2) "Fast Line-Profile Variability in the Spectra of O Stars" by Kholtygin, A. F., Monin, D. N., Surkov A. E., Fabrika, S. N. (2003)

        '''
        # loop over list of spectra
            # loop over wavelength 


class IueBinSys(Iue):
    '''
        This class will handle data from binary systems by adding an extra column containing the orbital phase of the observation.

        Usage:
            binSysName = IueBinSys() # new instance of the class
            binSysName.read_iue_data(fileList) # read data list
            binSysName.set_ephemeris(period, T0) # set ephemeris for the system
            binSysName.my_info.head() # show the resulting dataframe
    '''

    def __init__(self):
        pass

    def set_ephemeris(self, period, T0):
        '''
            Set the period and zero point for orbital phase calculation.
        '''

        if (self._my_list_ok == False):
            self._check_my_list()

        import astropy.units as u
        from astropy.time import Time, TimeDelta

        self.per = TimeDelta(period * u.day)
        self.T0 = Time(T0, format='jd')
        # ephemeris
        # T - T0
        deltaT = self.my_info['jdMid'] - self.T0.jd
        self.my_info['deltaT'] = deltaT
        # phase = mod( (T-T0) / per )           V  (T-T0) > 0
        # phase = mod(mod( (T-T0) / per ) + 1)  V  (T-T0) < 0
        self.my_info['phase'] = [ ( (dT/self.per.value)%1 if (dT > self.T0.jd) else ((dT/self.per.value)%1+1)%1 ) for dT in deltaT]

    def find_by_phase(self, start, end):
        '''
            Search dataset for data within a given interval in phase.

            Returns a list with a pandas dataframe for each element.
        '''
        import pandas as pd
        import numpy as np

        if (self._my_list_ok == False):
            self._check_my_list()

        s = start if start < end else end
        e = end if start < end else start

        # create a True/False mask
        mask = [(self.my_info.phase[i] >= s) & (self.my_info.phase[i] <= e) for i in range(len(self.my_info.DateObs))]

        # create list with indexes where mask == True
        mask_index = np.where(np.asarray(mask) == True)[0].tolist()

        # create series to be concatenated
        norders = [self.my_spec[mask_index[i]].shape[0] for i in range(len(mask_index))] # NORDERS IS VARIABLE!!
        nelements = len([item for item in self.my_info.jdObs.iloc[mask].values])
        # create a list with a pandas Series for each element
        s1 = [pd.Series(
            [self.my_info.DateObs.iloc[mask].values[i]] * norders[i],
            name='DateObs') for i in range(nelements)]
        s2 = [pd.Series(
            [self.my_info.aperture.iloc[mask].values[i]] * norders[i],
            name='aperture') for i in range(nelements)]
        s3 = [pd.Series(
            [self.my_info.jdObs.iloc[mask].values[i]] * norders[i],
            name='jdObs') for i in range(nelements)]
        s4 = [pd.Series(
            [self.my_info.jdMid.iloc[mask].values[i]] * norders[i],
            name='jdMid') for i in range(nelements)]
        s5 = [pd.Series(
            [self.my_info.exptime.iloc[mask].values[i]] * norders[i],
            name='exptime') for i in range(nelements)]

        # create a list with one valid data frame per element
        temp = [pd.concat([self.my_spec[mask_index[i]], s1[i], s2[i], s3[i], s4[i], s5[i]], axis=1) for i in range(len(mask_index))]

        # test for the existence of the key 'phase'
        if ('phase' in self.my_info):
            phase = [pd.Series(
                [self.my_info.phase.iloc[mask].values[i]] * norders[i], name='phase') for i in range(nelements)]
            # if present, concatenate it to the dataframe
            temp2 = [pd.concat([temp[i], phase[i]], axis=1) for i in range(len(temp))]
        else:
            # ignore it
            temp2 = temp

        # remove empty data
        res = []
        for item in temp2:
            if item.shape[0] != 0:
                res.append(item)

        # return a list with data that encompass the provided wavelength
        return res


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

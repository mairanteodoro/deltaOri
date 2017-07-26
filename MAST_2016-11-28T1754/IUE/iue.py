# !/usr/bin/env python3


class Iue(object):  # <- providing 'object' to comply with new style class!

    # flag for the existence of list; class attribute
    _my_list_ok = False

    def __init__(self, *arg):
        # arg == my_list
        # N.B.: the Iue class can be instantiated with or without providing a
        # list of the files to be read in and stored
        import pandas as pd

        self.my_list = ()
        self.my_info = pd.DataFrame()
        self.my_spec = []

        try:
            if arg.__len__() == 1:
                # if arg has been provided then read data immediately
                self.read_iue_data(arg[0])
            else:
                # otherwise, return instance of class
                return
        except RuntimeError as exception:
            print("Iue ERROR: {0}".format(exception))

    def _check_my_list(self):
        """
            Checks for the existence of IUE data.
        """
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
        except SystemExit as exception:
            print(exception)
            sys.exit()

    def read_iue_data(self, my_list: list):
        """
            This method appends 2 attributes to self: self.my_info and self.my_spec. They correspond to the info about
            the observation/data and the full spectrum, respectively. Both are saved in pandas data frame format.

            - In order to access the information about the files:
                self.my_info.head()

            - In order to access the full spectrum of a given observation:
                - first find the index of the desired file:
                idx = self.my_info[self.my_info.filename.str.contains("swp42749")].index[0]
                - now retrieve the wavelength and flux arrays:
                self.my_spec[idx].wave and self.my_spec[idx].flux
        """

        import pandas as pd
        import astropy.io.fits as fits

        self.my_list = my_list,  # <- this is a tuple!

        data_to_be_stored = [
            # [filename, PrimaryHDU, BinTableHDU]
            my_list,
            # [ list(fits.open(x))[0] for x in my_list ],
            # [ list(fits.open(x))[1] for x in my_list ],
            # info about observation
            #
            # APERTURE
            [fits.open(x)[0].header.get('APERTURE') for x in my_list],
            # DATE + TIME
            [fits.open(x)[0].header.get('SDATEOBS' if fits.open(x)[0].header.get('SDATEOBS') else 'LDATEOBS') +
             'T' + fits.open(x)[0].header.get('STIMEOBS' if fits.open(x)[0].header.get('STIMEOBS') else 'LTIMEOBS')
             for x in my_list],
            # EXPTIME
            [fits.open(x)[0].header.get('SEXPTIME' if fits.open(x)[0].header.get('SEXPTIME') else 'LEXPTIME') for x in
             my_list],
            # JD
            [fits.open(x)[0].header.get('SJD-OBS' if fits.open(x)[0].header.get('SJD-OBS') else 'LJD-OBS') for x in
             my_list],
            # MID JD
            [fits.open(x)[0].header.get('SJD-MID' if fits.open(x)[0].header.get('SJD-MID') else 'LJD-MID') for x in
             my_list]
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
        self.my_spec = [self.get_spec(x) for x in my_list]

    def get_spec(self, filename: str):
        """
            This method loops over the orders of a given file.

            Returns a dataframe with columns 'order', 'wave', 'flux', and 'filename'.
        """

        if self._my_list_ok is False:
            self._check_my_list()

        from astropy.table import Table
        import pandas as pd
        import numpy as np

        # find the filename's index
        file_index = self.find_index_of(filename)
        # read and return data in tabular format
        tb = Table.read(self.my_info.filename[file_index], hdu=1)
        wave, ripple, abs_cal, wave_min, wave_max, nu, noise, net = [], [], [], [], [], [], [], []

        order = [tb['ORDER'][i] for i in range(tb['ORDER'].shape[0])]
        for index in range(len(order)):
            # wavelength
            wave.append(tb['WAVELENGTH'][index] + (np.array(range((tb['RIPPLE'][index]).shape[0])) -
                        tb['STARTPIX'][index]) *
                        tb['DELTAW'][index])
            # wavelength intervals
            wave_min.append(wave[-1][0])
            wave_max.append(wave[-1][-1])
            # flux
            ripple.append(tb['RIPPLE'][index])
            abs_cal.append(tb['ABS_CAL'][index])
            nu.append(tb['QUALITY'][index])
            # this is **NOT** the flux-calibrated noise!
            # NOISE(flux-cal) = NOISE * (ABS_CAL / NET)
            noise.append(tb['NOISE'][index])
            net.append(tb['NET'][index])

        """
        Return a dataframe with the spectrum for each order:
             order   wave    fluxes    filename
             -----   -----   -----   -----
        0    #       []      []      string
        1    #       []      []      string
        """
        return pd.DataFrame({'filename': filename,
                             'order': order,
                             'wave': wave,
                             'ripple': ripple,
                             'abs_cal': abs_cal,
                             'waveMin': wave_min,
                             'waveMax': wave_max,
                             'quality': nu,
                             'noise': noise,  # <- this is not the calibrated noise!
                             'net_flux': net})

    def get_all_spec(self):
        """
            Return an array where each element is a pandas DataFrame containing data about 'order', 'wave', and 'flux'
            for the corresponding 'filename'.

            e.g.:
                iue = Iue()
                iue.read_iue_data(file)
                allSpec = iue.get_all_spec()
                allSpec[0].head()
                allSpec[0].filename
                allSpec[0].order, allSpec[0].wave, allSpec[0].flux
        """

        if not self._my_list_ok:
            self._check_my_list()

        # return an array where each element is a pandas DataFrame
        # containing 'filename', 'order', 'wave', and 'flux'
        return [self.get_spec(filename) for filename in self.my_list[0]]

    @staticmethod
    def _splice(w1, w2, f1, f2, q1, q2, n1=[], n2=[], net1=[], net2=[]):
        """
            Splice two spectra together taking into account data quality flag.

            If noise is provided (via n1, n2, net1, and net2), the returned value is flux-calibrated
            (i.e.: returned_noise(cgs) = noise(raw) * abs_cal(cgs) / net_flux(raw)).
            This could be useful for Temporal Variance analysis (default is NO NOISE).

            Adapted from:
            https://archive.stsci.edu/pub/iue/software/iuedac/procedures/mergeav.pro

            *** N.B.: w1[0] MUST BE < w2[0] ***

            E.g.: w1 = iue.my_spec[0].wave.iloc[0]
                  w2 = iue.my_spec[0].wave.iloc[1]
                  f1 = iue.my_spec[0].abs_cal.iloc[0]
                  f2 = iue.my_spec[0].abs_cal.iloc[1]
                  q1 = iue.my_spec[0].quality.iloc[0]
                  q2 = iue.my_spec[0].quality.iloc[1]
                ### OPTIONAL:
                  n1 = iue.my_spec[0].noise.iloc[0]
                  n2 = iue.my_spec[0].noise.iloc[1]
                  net1 = iue.my_spec[0].net_flux.iloc[0]
                  net2 = iue.my_spec[0].net_flux.iloc[1]
                ###
                  merged_spec = iue._splice(w1, w2, f1, f2, q1, q2, n1=n1, n2=n2, net1=net1, net2=net2)
        """
        from scipy import interpolate
        import numpy as np
        import astropy.units as u

        # units of the returned values
        # noinspection PyUnresolvedReferences
        wave_unit = u.Angstrom
        # noinspection PyUnresolvedReferences,PyUnresolvedReferences,PyUnresolvedReferences,PyUnresolvedReferences
        flux_unit = u.erg / u.s / (u.cm * u.cm) / u.Angstrom

        # considering only good quality flags
        good1, good2 = np.where(q1 >= 0)[0], np.where(q2 >= 0)[0]
        w1, f1, q1 = w1[good1], f1[good1], q1[good1]
        w2, f2, q2 = w2[good2], f2[good2], q2[good2]
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            n1, n2 = n1[good1], n2[good2]
            net1, net2 = net1[good1], net2[good2]

        # OVERLAPPING REGION
        # finding the overlapping region
        w1_overlap, f1_overlap, q1_overlap = w1[w1 >= w2[0]], f1[w1 >= w2[0]], q1[w1 >= w2[0]]
        w2_overlap, f2_overlap, q2_overlap = w2[w2 <= w1[-1]], f2[w2 <= w1[-1]], q2[w2 <= w1[-1]]
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            n1_overlap = n1[w1 >= w2[0]]
            n2_overlap = n2[w2 <= w1[-1]]
            net1_overlap = net1[w1 >= w2[0]]
            net2_overlap = net2[w2 <= w1[-1]]

        # NON OVERLAPPING REGION
        # indexes of the non-overlapping region
        non_overlap_indx = np.where(w1 >= w2[0])[0][0]
        # wave at the non-overlapping region
        w1_non_overlap, f1_non_overlap, q1_non_overlap = \
            w1[:non_overlap_indx], f1[:non_overlap_indx], q1[:non_overlap_indx]
        # w2_non_overlap, f2_non_overlap, q2_non_overlap = \
        #     w2[non_overlap_indx:], f2[non_overlap_indx:], q2[non_overlap_indx:]
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            n1_non_overlap = n1[:non_overlap_indx]
            # n2_non_overlap = n2[:non_overlap_indx]
            net1_non_overlap = net1[:non_overlap_indx]
            # net2_non_overlap = net2[:non_overlap_indx]

        # INTERPOLATION
        # determine the representation of (w2_overlap, f2_overlap)
        # using a cubic spline without smoothing (s=0)
        tck = interpolate.splrep(w2_overlap, f2_overlap, s=0)
        # interpolate the values of f2_overlap using the knots of w1_overlap
        f2new_overlap = interpolate.splev(w1_overlap, tck, der=0)
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            tck = interpolate.splrep(w2_overlap, n2_overlap, s=0)
            n2new_overlap = interpolate.splev(w1_overlap, tck, der=0)
            tck = interpolate.splrep(w2_overlap, net2_overlap, s=0)
            net2new_overlap = interpolate.splev(w1_overlap, tck, der=0)

        # SPLICING
        # take the average of f1_overlap and f2new_overlap
        f_overlap = np.average([f1_overlap, f2new_overlap], axis=0)
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            n_overlap = np.average([n1_overlap, n2new_overlap], axis=0)
            net_overlap = np.average([net1_overlap, net2new_overlap], axis=0)

        # OUTPUT
        # find the wavelengths of w2 that can be used to fill in the wavelengths of w1 that have q1 < 0
        w2_good = w2[np.where(w2 >= w1[-1])[0]]
        f2_good = f2[np.where(w2 >= w1[-1])[0]]
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            n2_good = n2[np.where(w2 >= w1[-1])[0]]
            net2_good = net2[np.where(w2 >= w1[-1])[0]]
        # use w1_overlap because f2_overlap was interpolated onto the w1_overlap array
        w3 = np.concatenate([w1_non_overlap, w1_overlap, w2_good])
        f3 = np.concatenate([f1_non_overlap, f_overlap, f2_good])
        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            n3 = np.concatenate([n1_non_overlap, n_overlap, n2_good])
            net3 = np.concatenate([net1_non_overlap, net_overlap, net2_good])

        if len(n1) != 0 and len(n2) != 0 and len(net1) != 0 and len(net2) != 0:
            # return wavelength, flux, and noise in cgs
            # return [w3, f3, n3*f3/net3]
            return {'wave': w3 * wave_unit, 'flux': f3 * flux_unit, 'noise': (n3 * f3 / net3) * flux_unit}
        else:
            # return wavelength and flux
            return {'wave': w3 * wave_unit, 'flux': f3 * flux_unit}

    def find_index_of(self, pattern):
        """
            Search .my_info for the index corresponding to the name of the file that matches up the provided pattern.
            The returned index can be used with anything within the .my_info dataframe.

            e.g.: index = iue.find_index_of("05722")
                  iue.exptime[index] # -> find exptime
                  iue.get_spec(iue.my_info.filename[index]) # -> get spectrum
        """

        if not self._my_list_ok:
            self._check_my_list()
        try:
            return self.my_info[self.my_info.filename.str.contains(pattern)].index[0]
        except IndexError:
            print("find_index_of exception: no match found.")
            # return None (type(returnedValue) = None)
            pass

    def find_by_date(self, start, end, wavelength=False):
        """
            Search dataframe for data within a given interval in date format.
            Optionally, if a wavelength is provided, filter the results by wavelength, as well
            (i.e. returns data filtered by date AND wavelength).

            Returns a list with a pandas dataframe for each element.
        """
        import pandas as pd
        import numpy as np

        if not self._my_list_ok:
            self._check_my_list()

        s = start if start < end else end
        e = end if start < end else start

        # create a True/False mask
        mask = [(self.my_info.DateObs[i] >= pd.Timestamp(s)) & (self.my_info.DateObs[i] <= pd.Timestamp(e)) for i in
                range(len(self.my_info.DateObs))]

        # create list with indexes where mask == True
        mask_index = np.where(np.asarray(mask))[0].tolist()

        # create series to be concatenated
        norders = [self.my_spec[mask_index[i]].shape[0] for i in range(len(mask_index))]  # NORDERS IS VARIABLE!!
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
        temp = [pd.concat([self.my_spec[mask_index[i]], s1[i], s2[i], s3[i], s4[i], s5[i]], axis=1) for i in
                range(len(mask_index))]

        # test for the existence of the key 'phase'
        if 'phase' in self.my_info:
            phase = [pd.Series(
                [self.my_info.phase.iloc[mask].values[i]] * norders[i], name='phase') for i in range(nelements)]
            # if present, concatenate it to the dataframe
            temp2 = [pd.concat([temp[i], phase[i]], axis=1) for i in range(len(temp))]
        else:
            # ignore it
            temp2 = temp

        # remove empty data
        res = []
        # for item in temp2:
        #     if item.shape[0] != 0:
        #         res.append(item)
        [res.append(item) for item in temp2 if item.shape[0] != 0]

        # return a list with data that encompass the provided wavelength
        if wavelength and wavelength > 0:
            # dataset filtered by phase and wavelength
            new_res = [
                df.query('waveMin <= {0} and waveMax >= {0}'.format(wavelength)) for df in res]
            # drop empty dataframes
            return [item for item in new_res if item.shape[0]]
        else:
            # dataset filtered by date only
            return res

    def find_by_wave(self, wave: float):

        """
            Find data that contains the provided wavelength.

            Returns a list with a pandas dataframe for each element.
        """
        import pandas as pd

        if self._my_list_ok is False:
            self._check_my_list()

        # create a True/False mask
        mask = [(self.my_spec[i].waveMin < wave) & (self.my_spec[i].waveMax > wave) for i in range(len(self.my_spec))]

        # create a list with one valid data frame per element
        # (note that j is not used as index; it is used to create a pd.Series of repeated values)
        # noinspection PyUnusedLocal,PyUnusedLocal,PyUnusedLocal,PyUnusedLocal,PyUnusedLocal
        temp = [pd.concat([self.my_spec[i][mask[i]],
                           pd.Series([self.my_info.DateObs[i] for j in range(len(self.my_spec[i][mask[i]]))],
                                     name='DateObs', index=self.my_spec[i][mask[i]].index),
                           pd.Series([self.my_info.aperture[i] for j in range(len(self.my_spec[i][mask[i]]))],
                                     name='aperture', index=self.my_spec[i][mask[i]].index),
                           pd.Series([self.my_info.jdObs[i] for j in range(len(self.my_spec[i][mask[i]]))],
                                     name='jdObs', index=self.my_spec[i][mask[i]].index),
                           pd.Series([self.my_info.jdMid[i] for j in range(len(self.my_spec[i][mask[i]]))],
                                     name='jdMid', index=self.my_spec[i][mask[i]].index),
                           pd.Series([self.my_info.exptime[i] for j in range(len(self.my_spec[i][mask[i]]))],
                                     name='exptime', index=self.my_spec[i][mask[i]].index)], axis=1) for i in
                range(len(self.my_spec))]

        # test for the existence of the key 'phase'
        if 'phase' in self.my_info:
            # if present, concatenate it to the dataframe
            # noinspection PyUnusedLocal
            temp2 = [pd.concat([temp[i],
                                pd.Series([self.my_info.phase[i] for j in range(len(self.my_spec[i][mask[i]]))],
                                          name='phase', index=self.my_spec[i][mask[i]].index)], axis=1) for i in
                     range(len(self.my_spec))]
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

    def tvs(spec_list: list, dispersion=None, continuum=None, S_par=0):
        """
        Uses the smoothed temporal variance (tvs) spectrum to analize light
        phase variations (LPV).

        Reference:
            (1) "Smoothed Temporal Variance Spectrum: weak line profile
            variations and NRP diagnostics" by Kholtygin, A. F., Sudnik,
            N. P. (2016);
            (2) "Fast Line-Profile Variability in the Spectra of O Stars" by
            Kholtygin, A. F., Monin, D. N., Surkov A. E., Fabrika, S. N. (2003)

        """
        # (1) loop over wavelength
        # (2) loop over list of spectra applying tvs
        from scipy import interpolate, ndimage
        import numpy as np

        # total number of spectra
        N = len(spec_list)
        const = 1 / (N - 1)

        smoothed_y = [ndimage.filters.gaussian_filter(item.get('flux'), S_par)
                      for item in spec_list]

        for i, item in enumerate(spec_list):
            item['flux'] = smoothed_y[i]

        # # normalize by the provided "continuum" region
        # flux_cnorm = [item.get('flux') /
        # item.get('flux')[(item.get('wave').value >= continuum[0]) & (
        #     item.get('wave').value <= continuum[1])].mean() for item in
        # spec_list]
        # noise_cnorm = [item.get('noise') /
        # item.get('flux')[(item.get('wave').value >= continuum[0]) & (
        #     item.get('wave').value <= continuum[1])].mean() for item in
        # spec_list]
        # # replace observed flux by continuum-normalized flux
        # for i, item in enumerate(spec_list):
        #     item['flux'] = flux_cnorm[i]
        #     item['noise'] = noise_cnorm[i]

        """ N.B.: the following is similar for both the TVS (Fullerton et al.
        1996) and tvs (Kholtygin and Sudnik 2016); the only difference being
        whether the input spec_list contains smoothed data or not, which can be
        defined by the S_par keyword (S_par=0 corresponds to non-smoothed data)
        """
        # sigma_ic =: noise in the continuum region = (S/N)^-1 = N/S =:
        #  stdev/mean
        sigma_ic = [
            item.get('noise')[
                (item.get('wave').value >= continuum[0]) &
                (item.get('wave').value <= continuum[1])].mean() /
            item.get('flux')[
                (item.get('wave').value >= continuum[0]) &
                (item.get('wave').value <= continuum[1])].mean()
            for item in spec_list]
        # mean noise for the entire dataset (sigma_zero)
        sigma_0 = ((1. / N) * np.sum([np.power(item.value, -2)
                                      for item in sigma_ic])) ** (-0.5)
        # sigma_i(lambda) =: (S/N)^-1 per pixel
        sigma_i = [item.get('noise') / item.get('flux') for item in spec_list]
        # weight of each spectrum based on the continuum SNR (Fullerton et al.
        # 1996)
        # w_i = (sigma_0 / sigma_ic)**2
        print('sigma_0: {0}; sigma_ic: {1}'.format(sigma_0, sigma_ic))

        """
            From Kholtygin and Sudnik (2016, pg 1609):

            "Note that the variation of the parameter S and the ratio S/N
            extremely weakly influence the value of
            q(lambda, S). This result holds for different values of the line
            profile parameters and the ratio S/N.
            Therefore for all S we can put q(lambda, S) approximately equals to
            q(lambda, 0). Thus for all values
            of S it is possible to use the simple formula (23) instead of
            making the numerical integration in
            equation (25)."
        """
        q_i = (sigma_0 / sigma_i) ** 2
        # print('q_i[0]: {}'.format(q_i[0]))
        # residual matrix
        d_i = [item.get('flux') - item.get('flux').mean()
               for item in spec_list]
        # xi_i^2
        xi_i2 = [q_i[i] * np.power(d_i[i], 2) for i in range(len(q_i))]

        # TODO: perform the clipping/interpolation BEFORE any math op

        # clip the spectra to the same wavelength range
        minima = [np.min(item.get('wave').value) for item in spec_list]
        maxima = [np.max(item.get('wave').value) for item in spec_list]
        x_limits = (np.max(minima), np.min(maxima))
        print('minima and maxima: {0} and {1}'.format(minima, maxima))
        print('x_limits: {0} -- {1}'.format(x_limits[0], x_limits[1]))
        x_clipped = [
            item.get('wave')[(item.get('wave').value >= x_limits[0]) &
                             (item.get('wave').value <= x_limits[1])]
            for item in spec_list]
        y_clipped = [
            product[(spectrum.get('wave').value >= x_limits[0]) &
                    (spectrum.get('wave').value <= x_limits[1])]
            for product, spectrum in zip(xi_i2, spec_list)]
        # new x grid where the interpolation will take place
        # x_new_min, x_new_max = np.max(minima), np.min(maxima)
        # find function shape for each spectrum
        func_interp = [
            interpolate.splrep(x_clipped[i].value, y_clipped[i].value)
            for i in range(len(x_clipped))]
        # new wavelength grid with provided dispersion
        x_new = np.arange(x_limits[0], x_limits[1], dispersion)
        # apply the interpolation function onto the new wavelength grid
        y_new = [interpolate.splev(x_new, func_interp[i])
                 for i in range(len(func_interp))]
        # import ipdb; ipdb.set_trace()
        # gaussian smooth (maybe this should be done at the beginning?)
        # tvs = ndimage.filters.gaussian_filter(const *
        #       np.sum(np.power(y_new, 2), axis=0), S_par)
        tvs = const * np.sum(np.power(y_new, 2), axis=0)

        return {'wave': x_new, 'tvs': tvs}


class IueBinSys(Iue):
    """
        This class will handle data from binary systems by adding an extra
        column containing the orbital phase of the observation.

        Usage:
            binSysName = IueBinSys() # new instance of the class
            binSysName.read_iue_data(fileList) # read data list
            binSysName.set_ephemeris(period, T0) # set ephemeris for the system
            binSysName.my_info.head() # show the resulting dataframe
    """

    def __init__(self, *arg):
        # calling parent's __init__()
        super().__init__(*arg)
        self.per = 0
        self.T0 = 0

    def set_ephemeris(self, period, T0):
        """
            Set the period and zero point for orbital phase calculation.
        """

        if self._my_list_ok is False:
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
        self.my_info['phase'] = [
            ((dT / self.per.value) % 1 if (dT > self.T0.jd) else
             ((dT / self.per.value) % 1 + 1) % 1) for dT in deltaT]

    def find_by_phase(self, start, end, wavelength=False):
        """
            Search dataset for data within a given interval in phase.
            Optionally, if a wavelength is provided, filter the results by
            wavelength, as well (i.e. returns data filtered by phase AND
            wavelength).

            Returns a list with a pandas dataframe for each element.
        """
        import pandas as pd
        import numpy as np

        if self._my_list_ok is False:
            self._check_my_list()

        s = start if start < end else end
        e = end if start < end else start

        # create a True/False mask
        mask = [(self.my_info.phase[i] >= s) & (self.my_info.phase[i] <= e)
                for i in range(len(self.my_info.DateObs))]

        # create list of indexes where mask == True
        mask_index = np.where(np.asarray(mask))[0].tolist()

        # create series to be concatenated
        # CAREFUL HERE; NORDERS IS VARIABLE!!
        norders = [self.my_spec[mask_index[i]].shape[0]
                   for i in range(len(mask_index))]
        nelements = len([item
                         for item in self.my_info.jdObs.iloc[mask].values])
        # create a list of pandas Series for each element
        s1 = [pd.Series(
            # creating norders copies of DateObs value
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

        # create a list of one valid data frame per element
        temp = [pd.concat([self.my_spec[mask_index[i]],
                           s1[i], s2[i], s3[i], s4[i], s5[i]], axis=1)
                for i in range(len(mask_index))]

        # test for the existence of the key 'phase'
        if 'phase' in self.my_info:
            phase = [pd.Series(
                [self.my_info.phase.iloc[mask].values[i]] * norders[i],
                name='phase') for i in range(nelements)]
            # if present, concatenate it to the dataframe
            temp2 = [pd.concat([temp[i], phase[i]], axis=1)
                     for i in range(len(temp))]
        else:
            # ignore it
            temp2 = temp

        # remove empty data
        res = []
        # for item in temp2:
        #     if item.shape[0] != 0:
        #         res.append(item)
        [res.append(item) for item in temp2 if item.shape[0] != 0]

        # return a list with data obtained between the provided phases
        if wavelength and wavelength > 0:
            # dataset filtered by phase and wavelength
            new_res = [df.query('waveMin <= {0} and waveMax >= \
                                {0}'.format(wavelength)) for df in res]
            # drop empty dataframes
            return [item for item in new_res if item.shape[0]]
        else:
            # dataset filtered by phase only
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

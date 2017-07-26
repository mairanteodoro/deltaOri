# !/usr/bin/env python3

###
import matplotlib.pyplot as plt
import json as json
from iue import Iue, IueBinSys
import numpy as np

# parse filenames into a list format
fileList = json.load(open("./iueData.json"))


# TWO WAYS TO INSTANTIATE THE Iue CLASS:

# - 1:
# instantiate the Iue class for general IUE data
iue = Iue()
# and do the actual reading of the data and append it to the instance
iue.read_iue_data(fileList['filename'])

# - 2:
# instantiate the Iue class providing the fileList
iue2 = Iue(fileList['filename'])


### DATA MANIPULATION
# query data by start and end dates
sel1_by_date = iue.find_by_date('1981-01', '1983-12')
# it doesn't matter the order of the dates
sel2_by_date = iue.find_by_date('1983-12', '1981-01')

# get the spectrum from a single file
singleSpec1 = iue.get_spec("./DATA/lwr02241.mxhi")
# read and store all spectra
allSpec1 = iue.get_all_spec()

# extract the data (df format) from file lwr02992.mxhi
spec0 = allSpec1[iue.find_index_of('02992')]
# extract information (df format) about file lwr02992.mxhi
info0 = iue.my_info.loc[iue.find_index_of('02992')]

# query by wavelength
sel1_by_wave = iue.find_by_wave(1830)


# PLOTTING
# plot all the orders
fig, ax = plt.subplots()

[ax.plot(spec0.wave.iloc[i], spec0.abs_cal.iloc[i], marker='')
 for i in range(len(spec0.wave))]

fig.savefig('./RESULTS/spec.pdf')

plt.clf()
plt.close('all')


def setup():

    ###
    # delta Ori
    ###
    import json as json
    from iue import Iue, IueBinSys

    # parse filenames into a list format
    fileList = json.load(open("./iueData.json"))

    # instantiate the IueBinSys class for binary systems
    delOri = IueBinSys()
    delOri.read_iue_data(fileList['filename'])

    # setting the ephemeris for the binary system
    #
    # ephemeris for delta Ori from Corcoran et al. 2015
    per, T0 = 5.732436, 2456295.674
    # set the ephemeris for delta Ori
    delOri.set_ephemeris(per, T0)

    return delOri


def do_it(wavelength: float=0,
          continuum: list=[],
          phase1: float =0,
          phase2: float =1
          ):

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u

    wave_unit = u.Angstrom
    flux_unit = u.erg / u.s / (u.cm * u.cm) / u.Angstrom

    # query data by phase and wavelength
    # this will store **two adjacent orders** from all the
    # data within the provided phase interval
    selected_data = delOri.find_by_phase(phase1, phase2, wavelength=wavelength)

    # import ipdb; ipdb.set_trace()
    spec_list = []
    for i, item in enumerate(selected_data):
        # get the number of orders found
        if item.shape[0] == 1:
            # only one order returned
            spec_list.append({'wave': item.wave.values[0] * wave_unit,
                              'flux': item.abs_cal.values[0] * flux_unit,
                              'noise': (item.noise.values[0] *
                                        item.abs_cal.values[0] /
                                        item.net_flux.values[0]) * flux_unit})
        else:
            # select the orders to merge using selection by position
            # (.iloc[] can be used too)
            w1, w2 = zip(item.wave)
            f1, f2 = zip(item.abs_cal)
            q1, q2 = zip(item.quality)
            n1, n2 = zip(item.noise)  # this is used for the SNR estimate
            net1, net2 = zip(item.net_flux)  # this is used for SNR estimate
            # merge the two orders (the returned noise is in cgs units)
            spec_list.append(delOri._splice(w1[0], w2[0],
                                            f1[0], f2[0],
                                            q1[0], q2[0],
                                            n1[0], n2[0],
                                            net1[0], net2[0]
                                            ))
    disp = 0.03  # Angstrom
    # cont_region = continuum # continuum region for normalization
    # and SNR estimate
    S_parameter = [0, 1, 5, 10]  # S parameter

    # smoothed TVS
    smTVS = [delOri.tvs(
            spec_list=spec_list,
            dispersion=disp,
            continuum=continuum,
            S_par=S
            ) for S in S_parameter]

    #
    # *** DONT FORGET TO SET plt.ion() ***
    plt.close('all')
    plt.ion()
    #
    fig, ax = plt.subplots(2, 1, sharex=True)
    # plot spectra
    [ax[0].plot(spec_list[i]['wave'], spec_list[i]['flux'],
                label=r'epoch {}'.format(i))
     for i in range(len(spec_list))]
    ax[0].set_ylabel(r'flux [erg\,s$^{-1}$\,cm$^{-2}$\,\AA$^{-1}$]')
    ax[0].legend(loc='best')

    # plot tvs for different values of S
    ax[1].plot(smTVS[0]['wave'], (smTVS[0]['tvs']),
               label=r'$\sigma$ = {0:0.2f}~\AA'.format(S_parameter[0] * disp))
    ax[1].set_xlabel(r'Wavelength [\AA]')
    ax[1].set_ylabel(r'(TVS)$^{1/2}$')
    ax[1].legend(loc='best')

    # ax[1].axhline(y=1.4176, ls='dashed')

    import ipdb; ipdb.set_trace()


delOri = setup()
do_it(wavelength=1400, phase1=0., phase2=0.2, continuum=[1385, 1387])
do_it(wavelength=1550, phase1=0., phase2=0.2, continuum=[1542, 1544])



### query data by date interval (the start/end date order does not matter)
sel_1 = delOri.find_by_date('1981-01', '1983-12')

### query data by phase interval (the start/end phase order does not matter)
sel_2 = delOri.find_by_phase(0.45, 0.55)

### query by wavelength only
sel_NV = delOri.find_by_wave(1240)
sel_CIV = delOri.find_by_wave(1550)
sel_SiIV = delOri.find_by_wave(1400)
# once we have selected data by wavelength, we can search it by phase using a two-step procedure
phase_interval = [0.45, 0.55]
# firstly, select all data that complies with the wavelength constraint
selected_data = [df.query('phase >= {0} and phase <= {1}'.format(phase_interval[0], phase_interval[1]))
for df in sel_NV]
sel_CIV_phase = [df.query('phase >= {0} and phase <= {1}'.format(phase_interval[0], phase_interval[1]))
for df in sel_CIV]
sel_SiIV_phase = [df.query('phase >= {0} and phase <= {1}'.format(phase_interval[0], phase_interval[1]))
for df in sel_SiIV]
# secondly, drop all dataframe that returned empty from the previous line
selected_data = [item for item in selected_data if item.shape[0]]
sel_CIV_phase = [item for item in sel_CIV_phase if item.shape[0]]
sel_SiIV_phase = [item for item in sel_SiIV_phase if item.shape[0]]
# to prove that both approaches (query by phase and then wavelength or by wavelength and then phase)
# return the same results, we can compare them:
# (will return [True, True, True, True])
[selected_data[i].equals(selected3_NV[i]) for i in range(len(selected_data))]

# plot the distribution of phase coverage
# plt.ion()
# delOri.myInfo.hist(column='phase')

# # plot all the spectra and orders
# fig, ax = plt.subplots()
#
# data = sel_SiIV
# colors = ['red', 'green']
#
# [[ax.plot(data[j].wave.iloc[i], data[j].abs_cal.iloc[i], marker='', color=colors[i], lw=0.1)
# for i in range(len(data[j].wave))] for j in range(len(data))]
#
# ax.set_xlabel(u'Wavelength [\AA]')
# ax.set_ylabel(u'abs\_cal [erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')
#
# fig.savefig('./RESULTS/spec_SiIV.svg')
#
# plt.clf()
# plt.close('all')

























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

import matplotlib.pyplot as plt
from scipy import interpolate
# import statistics as stats
# import sys
# this backend allows for multiple pages in a PDF file
from matplotlib.backends.backend_pdf import PdfPages

# loop over the files
for index in list(range(iue.my_spec.__len__())):

    currSpec = iue.my_spec[index]

    # create a single PDF file with multiple pages
    pdf_pages = PdfPages(
        './RESULTS/{0}_orders.pdf'.format(
            (set(currSpec.filename).pop()).split('/')[-1]
        )
    )

    # loop over the orders for each file
    for i in list(range(currSpec.order.__len__())):
        print('{0}:  i = {1}'.format(
                (set(currSpec.filename).pop()).split('/')[-1],
                i,
            )
        )
        try:
            w1, w2 = currSpec.wave[0+i:2+i].values
            r1, r2 = currSpec.ripple[0+i:2+i].values
            a1, a2 = currSpec.abs_cal[0+i:2+i].values

            # stats
            # r1_mean, r1_sdev = stats.mean(r1), stats.stdev(r1)
            # a1_mean, a1_sdev = stats.mean(a1), stats.stdev(a1)
            # r2_mean, r2_sdev = stats.mean(r2), stats.stdev(r2)
            # a2_mean, a2_sdev = stats.mean(a2), stats.stdev(a2)

            subset_w1, subset_r1, subset_a1 = w1[w1 >= w2[0]], r1[w1 >= w2[0]], a1[w1 >= w2[0]]
            subset_w2, subset_r2, subset_a2 = w2[w2 <= w1[-1]], r2[w2 <= w1[-1]], a2[w2 <= w1[-1]]

            tck_r1 = interpolate.splrep(subset_w1, subset_r1, s=0)
            subset_r1new = interpolate.splev(subset_w2, tck_r1, der=0)
            tck_a1 = interpolate.splrep(subset_w1, subset_a1, s=0)
            subset_a1new = interpolate.splev(subset_w2, tck_a1, der=0)

            # create a new page
            fig, ax = plt.subplots(
                nrows=2, ncols=1,
                sharex=True, squeeze=True,
                figsize=(8.27, 11.69), dpi=100
            )

            plt.subplots_adjust(hspace=0.05)

            ax[1].set_xlabel(r'Wavelength (\AA)')

            ax[0].set_title('{0}: orders {1} and {2}'.format(
                    (set(currSpec.filename).pop()).split('/')[-1],
                    currSpec.order[0+i:2+i].values[0],
                    currSpec.order[0+i:2+i].values[1],
                )
            )

            ax[0].annotate('RIPPLE',
                xy=(0.05,0.95),
                xycoords='axes fraction',
                size=14, ha='left', va='top',
            )

            ax[1].annotate('ABS\_CAL',
                xy=(0.05,0.95),
                xycoords='axes fraction',
                size=14, ha='left', va='top',
            )

            ax[0].plot(w1, r1, marker='')
            ax[0].plot(w2, r2, marker='', color='red')
            ax[0].scatter(subset_w2, subset_r1new, s=5, facecolors='none', edgecolors='k', zorder=10)

            ax[1].plot(w1, a1, marker='')
            ax[1].plot(w2, a2, marker='', color='red')
            ax[1].scatter(subset_w2, subset_a1new, s=5, facecolors='none', edgecolors='k', zorder=10)

            pdf_pages.savefig(fig)

        except:
            pass

    pdf_pages.close()
    plt.close('all')

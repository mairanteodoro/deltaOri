# !/usr/bin/env python3

###
import matplotlib.pyplot as plt
import json as json
from iue import Iue, IueBinSys
import numpy as np

# parse filenames into a list format
fileList = json.load( open( "./iueData.json" ) )


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

[ax.plot(spec0.wave.iloc[i], spec0.abs_cal.iloc[i], marker='') for i in range(len(spec0.wave))]

fig.savefig('./RESULTS/spec.pdf')

plt.clf()
plt.close('all')



###
# delta Ori
###
# instantiate the IueBinSys class for binary systems
delOri = IueBinSys()
delOri.read_iue_data(fileList['filename'])

# ephemeris for delta Ori from Corcoran et al. 2015
per, T0 = 5.732436, 2456295.674
# set the ephemeris for delta Ori
delOri.set_ephemeris(per, T0)


# query data by start and end dates
selected1 = delOri.find_by_date('1981-01', '1983-12')
# or
selected2 = delOri.find_by_date('1983-12', '1981-01')

# query data by start and end phase
selected3 = delOri.find_by_phase(0.1, 0.2)
# or
selected4 = delOri.find_by_phase(0.2, 0.1)

# query by wavelength
sel_NV = delOri.find_by_wave(1240)
sel_CIV = delOri.find_by_wave(1550)
sel_SiIV = delOri.find_by_wave(1400)

# select orders to merge
w1 = sel_NV[0].wave.iloc[0]
w2 = sel_NV[0].wave.iloc[1]
f1 = sel_NV[0].abs_cal.iloc[0]
f2 = sel_NV[0].abs_cal.iloc[1]
q1 = sel_NV[0].quality.iloc[0]
q2 = sel_NV[0].quality.iloc[1]
# merge the two orders
merged = delOri._splice(w1, w2, f1, f2, q1, q2)

# plot the distribution of phase coverage
# plt.ion()
# delOri.myInfo.hist(column='phase')

# plot all the spectra and orders
fig, ax = plt.subplots()

data = sel_SiIV
colors = ['red', 'green']

[[ax.plot(data[j].wave.iloc[i], data[j].abs_cal.iloc[i], marker='', color=colors[i], lw=0.1) for i in range(len(data[j].wave))] for j in range(len(data))]

ax.set_xlabel(u'Wavelength [\AA]')
ax.set_ylabel(u'abs\_cal [erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')

fig.savefig('./RESULTS/spec_SiIV.svg')

plt.clf()
plt.close('all')

























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
import statistics as stats
import sys
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

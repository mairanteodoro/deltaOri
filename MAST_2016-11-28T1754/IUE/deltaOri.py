# !/usr/bin/env python

###
import matplotlib.pyplot as plt
import json as json
from iue import Iue, IueBinSys

# parse filenames into a list format
fileList = json.load( open( "./iueData.json" ) )



# instantiate the Iue class for general IUE data
iue = Iue()
# do the actual reading of the data and append it to the instance
iue.readIueData(fileList['filename'])

### DATA MANIPULATION
# query data by start and end dates
sel1 = iue.findByDate('1981-01', '1983-12')
# it doesn't matter the order of the dates
sel2 = iue.findByDate('1983-12', '1981-01')

# get the spectrum from a single file
singleSpec1 = iue.getSpec("./DATA/lwr02241.mxhi")
# read and store all spectra
allSpec1 = iue.getAllSpec()

# extract the data (df format) from file lwr02992.mxhi
spec0 = allSpec1[iue.findIndexOf('02992')]
# extract information (df format) about file lwr02992.mxhi
info0 = iue.myInfo.loc[iue.findIndexOf('02992')]

# plot all the orders
fig, ax = plt.subplots()

[ax.plot(spec0.wave[i], spec0.flux[i], marker='') for i in range(len(spec0.wave))]

fig.savefig('./RESULTS/spec.pdf')

plt.clf()
plt.close('all')




# instantiate the IueBinSys class for binary systems
delOri = IueBinSys()
delOri.readIueData(fileList['filename'])

# ephemeris for delta Ori from Corcoran et al. 2015
per, T0 = 5.732436, 2456295.674
# set the ephemeris for delta Ori
delOri.setEphemeris(per, T0)

# query data by start and end dates
selected1 = delOri.findByDate('1981-01', '1983-12')
# or
selected2 = delOri.findByDate('1983-12', '1981-01')

# query data by start and end phase
selected3 = delOri.findByPhase(0.1, 0.2)
# or
selected4 = delOri.findByPhase(0.2, 0.1)

# plot the distribution of phase coverage
plt.ion()
delOri.myData.hist(column='phase')


















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

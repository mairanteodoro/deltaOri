###
import json as json
from iue import Iue

# ephemeris for delta Ori from Corcoran et al. 2015
per, T0 = 5.732436, 2456295.674

# parse filenames into a list format
fileList = json.load( open( "./iueData.json" ) )

# instantiate the Iue class and set ephemeris
delOri = Iue(per, T0)

# do the actual reading of the data and append it to the instance
delOri.readIueData( fileList['filename'] )

### MANIPULATION OF DATA
# query data by start and end dates
selected1 = delOri.findByDate('1981-01', '1983-12')
# it doesn't matter the order of the dates
selected2 = delOri.findByDate('1983-12', '1981-01')

# query data by start and end phase
selected3 = delOri.findByPhase(0.1, 0.2)
# or
selected4 = delOri.findByPhase(0.2, 0.1)












# # read and store all spectra
# # using a JSON file
# allSpec = res.getAllSpec("./iueData.json")
# # using the output from json.load()
# allSpec2 = res.getAllSpec(fileList['filename'])




# get the distribution of phase coverage
# plt.ion()
# res.myData.hist(column='phase')




# get all the orders for a file
# spec = res.createSpec("./DATA/lwr02241.mxhi")

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

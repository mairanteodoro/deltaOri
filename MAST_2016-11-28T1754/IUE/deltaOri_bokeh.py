# !/usr/bin/env python3

###
import matplotlib.pyplot as plt
import json as json
from iue import Iue, IueBinSys
from bokeh.plotting import output_file, show, figure, curdoc
from bokeh.client import push_session
from bokeh.mpl import to_bokeh
import numpy as np

# parse filenames into a list format
fileList = json.load( open( "./iueData.json" ) )

# instantiate the IueBinSys class for binary systems
delOri = IueBinSys()
delOri.read_iue_data(fileList['filename'])

# ephemeris for delta Ori from Corcoran et al. 2015
per, T0 = 5.732436, 2456295.674
# set the ephemeris for delta Ori
delOri.set_ephemeris(per, T0)

test = delOri.find_by_date('1981-02-15', '1981-02-20')

sel_by_date1 = delOri.find_by_date('1978-04', '1978-05')
sel_by_phase1 = delOri.find_by_phase(0.1, 0.2)
sel_by_wave1 = delOri.find_by_wave(1240)


d1, d2 = sel_by_wave1[0].iloc[0], sel_by_wave1[0].iloc[1]




# plot all the spectra and orders
fig, ax = plt.subplots()

colors = ['red', 'green']

[[ax.plot(sel_by_wave1[j].wave.iloc[i], sel_by_wave1[j].abs_cal.iloc[i], marker='', color=colors[i], lw=0.1) for i in range(len(sel_by_wave1[j].wave))] for j in range(len(sel_by_wave1))]

ax.set_xlabel(u'Wavelength [\AA]')
ax.set_ylabel(u'abs\_cal [erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')




# USING BOKEH
TOOLS = 'pan, box_zoom, wheel_zoom, hover, resize, undo, redo, reset, save'
# x = sel_by_wave1[0].wave.iloc[0]
# y = sel_by_wave1[0].abs_cal.iloc[0]
# # output_file('delOri.html')
# p = figure(tools=TOOLS)
#
# [p.line(x, y, color=color) for x, y, color in zip(sel_by_wave1[0].wave, sel_by_wave1[0].abs_cal, ['red', 'green'])]
#
# session = push_session(curdoc())
#
# session.show(p)
# session.loop_until_closed()
output_file('delOri.html')
show(to_bokeh(fig, tools=TOOLS))

fig.savefig('./RESULTS/spec_l1240.pdf')

plt.clf()
plt.close('all')

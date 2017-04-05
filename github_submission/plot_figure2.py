import matplotlib.pylab as pt
import numpy
import matplotlib.cm as cm

import scipy.stats
from matplotlib import colors

import matplotlib

matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.size'] = 10.0
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['lines.markersize'] = 3.5
matplotlib.rcParams['lines.linewidth'] = .7
matplotlib.rcParams['legend.fontsize'] = 8.0
matplotlib.rcParams['axes.linewidth']= .5
matplotlib.rcParams['patch.linewidth']= .5

##This script plots the heritability statistics on initial fitness and adaptability with 2-sigma confidence intervals, as calculated
##and recorded in 'write_heritability_statistics_table.py'.

filename_in = 'data/heritability_statistics_with_differences_2_15_2017.csv'
file = open(filename_in,'r')

stats_dict = {}
for line in file:
	if line.startswith('#'):
		stat_name = line.strip()[1:]
		print stat_name
		if '-' not in stat_name:
			stats_dict[stat_name] = {}
	elif '-' not in stat_name:
		line_list = line.strip().split(',')
		stats_dict[stat_name][line_list[0]] = [float(l) for l in line_list[1:]] #mean, lower, and upper bound for statistic
		stats_dict[stat_name][line_list[0]][1] = stats_dict[stat_name][line_list[0]][0] - stats_dict[stat_name][line_list[0]][1] #Write interval rather than the lower bound, because this is the format for matplotlib's error bars
		stats_dict[stat_name][line_list[0]][2] = stats_dict[stat_name][line_list[0]][2] - stats_dict[stat_name][line_list[0]][0]

file.close()

##Colorbrewer2 pastel set 2

color1b = numpy.array((166,206,227),dtype='float')/255.
color2b = numpy.array((31,120,180),dtype='float')/255.
color3b = numpy.array((178,223,138),dtype='float')/255.
color4b = numpy.array((51,160,44),dtype='float')/255.
color5b = numpy.array((251,154,153),dtype='float')/255.

###Colorbrewer2 set 3

color1c = numpy.array((27,158,119),dtype='float')/255.
color2c = numpy.array((217,95,2),dtype='float')/255.
color3c = numpy.array((117,112,179),dtype='float')/255.
color4c = numpy.array((231,41,138),dtype='float')/255.
color5c = numpy.array((102,166,30),dtype='float')/255.

##Colorbrewer set 4
color1d = numpy.array((141,211,199),dtype='float')/255.
color2d = numpy.array((255,255,179),dtype='float')/255.
color3d = numpy.array((190,186,218),dtype='float')/255.
color4d = numpy.array((251,128,114),dtype='float')/255.
color5d = numpy.array((128,177,211),dtype='float')/255.


##Colorbrewer set 5

color1e = numpy.array((166,97,26),dtype='float')/255.
color2e = numpy.array((223,194,125),dtype='float')/255.
color3e = numpy.array((245,245,245),dtype='float')/255.
color4e = numpy.array((128,205,193),dtype='float')/255.
color5e = numpy.array((1,133,113),dtype='float')/255.


fig, (ax1, ax2) = pt.subplots(1,2, figsize = (8,4))
b1 = ax1.bar([0,4],[stats_dict['H^2']['Initial fitness YPD at 30C'][0], stats_dict['H^2']['Initial fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['H^2']['Initial fitness YPD at 30C'][1:], stats_dict['H^2']['Initial fitness SC at 37C'][1:]]).T, color=color1e, ecolor='k',capsize=0)
b2 = ax1.bar([1,5],[stats_dict['h^2']['Initial fitness YPD at 30C'][0], stats_dict['h^2']['Initial fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['h^2']['Initial fitness YPD at 30C'][1:], stats_dict['h^2']['Initial fitness SC at 37C'][1:]]).T, color=color2e, ecolor='k',capsize=0)
b3 = ax1.bar([2,6],[stats_dict['r^2, QTL models']['Initial fitness YPD at 30C'][0], stats_dict['r^2, QTL models']['Initial fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['r^2, QTL models']['Initial fitness YPD at 30C'][1:], stats_dict['r^2, QTL models']['Initial fitness SC at 37C'][1:]]).T, color=color3e, ecolor='k',capsize=0)
ax1.legend([b1, b2, b3], ['Broad-sense heritability','Narrow-sense heritability','QTL model, var. explained'])
ax1.set_ylim(0,1.3)
ax1.set_xlim(-.2,7)
ax1.set_xticks([1.5, 5.5])
ax1.set_xticklabels(['Fitness at OT','Fitness at HT'])
b1 = ax2.bar([0,6],[stats_dict['H^2']['Delta fitness YPD at 30C'][0], stats_dict['H^2']['Delta fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['H^2']['Delta fitness YPD at 30C'][1:], stats_dict['H^2']['Delta fitness SC at 37C'][1:]]).T, color=color1e, ecolor='k',capsize=0)
b2 = ax2.bar([1,7],[stats_dict['h^2']['Delta fitness YPD at 30C'][0], stats_dict['h^2']['Delta fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['h^2']['Delta fitness YPD at 30C'][1:], stats_dict['h^2']['Delta fitness SC at 37C'][1:]]).T, color=color2e, ecolor='k',capsize=0)
b3 = ax2.bar([2,8],[stats_dict['r^2, QTL models']['Delta fitness YPD at 30C'][0], stats_dict['r^2, QTL models']['Delta fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['r^2, QTL models']['Delta fitness YPD at 30C'][1:], stats_dict['r^2, QTL models']['Delta fitness SC at 37C'][1:]]).T, color=color3e, ecolor='k',capsize=0)
b4 = ax2.bar([3,9],[stats_dict['r^2, fitness models']['Delta fitness YPD at 30C'][0], stats_dict['r^2, fitness models']['Delta fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['r^2, fitness models']['Delta fitness YPD at 30C'][1:], stats_dict['r^2, fitness models']['Delta fitness SC at 37C'][1:]]).T, color=color4e, ecolor='k',capsize=0)
b5 = ax2.bar([4,10],[stats_dict['r^2, combined models']['Delta fitness YPD at 30C'][0], stats_dict['r^2, combined models']['Delta fitness SC at 37C'][0]], yerr = numpy.array([stats_dict['r^2, combined models']['Delta fitness YPD at 30C'][1:], stats_dict['r^2, combined models']['Delta fitness SC at 37C'][1:]]).T, color=color5e, ecolor='k',capsize=0)
ax2.legend([b1, b2, b3, b4, b5], ['Broad-sense heritability','Narrow-sense heritability','QTL model','Fitness model','Combined model'],ncol=1)
ax2.set_xticks([2.5, 8.5])
ax2.set_ylim(0,1)
ax2.set_xticklabels(['Fitness gains at OT','Fitness gains at HT'])
ax2.set_xlim(-.2, 11)

pt.savefig('bar_chart_2_23_2017.pdf',bbox_inches='tight')
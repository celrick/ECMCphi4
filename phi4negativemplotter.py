import matplotlib.pyplot as plt
import numpy as np
xs = [6,8,10,12,14]

"""
#av action
#fill in these with data
mcmcs = [0.2476, 0.2512, 0.2669,0.3252,0.3354]
ecmcs = [0.2453,0.2326,0.2376,0.2307,0.2336]
mcmcerrors = [0.0156,0.0228,0.0401,0.1002,0.1107]
ecmcerrors = [0.0141,0.0051,0.0105,0.0048,0.0085]
observable = 'Average Action per site'
"""
"""
#obs1
mcmcs = [0.2312,0.2474,0.2566,0.2629,0.2669]
ecmcs =[0.2305,0.2465,0.2567,0.2626,0.2662]
mcmcerrors = [6.6914E-4,6.7633E-4,7.1924E-4,0.0013,0.0012]
ecmcerrors = [5.9365E-4,4.5277E-4,4.8374E-4,3.5537E-4,4.6088E-4]
observable = 'Observable 1'
"""
"""
#obs 2
mcmcs = [0.0668,0.0724,0.0759,0.0783,0.0797]
ecmcs = [0.0668,0.0723,0.0758,0.0780,0.0794]
mcmcerrors = [1.5502E-4,1.7319E-4,2.0608E-4,3.7829E-4,3.5526E-4]
ecmcerrors = [1.5286E-4,1.2234E-4,1.2715E-4,9.1591E-5,1.1836E-4]
observable = 'Observable 2'
"""
"""
#obs 3
mcmcs = [0.0035,0.0022,0.0015,0.0012,8.547E-4]
ecmcs = [0.0035,0.0023,0.0015,0.0011,8.5493E-4]
mcmcerrors = [4.8401E-5,3.0756E-5,2.1988E-5,1.6933E-5,1.2206E-5]
ecmcerrors = [4.7795E-5,3.2409E-5,2.144E-5,1.5414E-5,1.1933E-5]
observable = 'Observable 3'
"""
"""
#ac action
mcmcs = [0.7312,0.7209,1.2473,3.2737,4.6607]
ecmcs = [0.5021,0.5001,0.4993,0.5011,0.5005]
mcmcerrors = [0.0319,0.0315,0.0759,0.2945,0.4829]
ecmcerrors = [0.0142,0.0100,0.01,0.0142,0.0142]
observable = 'Integrated Autocorrelation time of Average Action'
"""

#Ac obs 1
mcmcs = [0.5546,0.5861,0.8129,1.5435,2.0886]
ecmcs = [0.5144,0.4901,0.4884,0.5032,0.4974]
mcmcerrors = [0.0190,0.0232,0.0388,0.1022,0.1530]
ecmcerrors = [0.0178,0.0098,0.0098,0.0142,0.0141]
observable = 'Integrated Autocorrelation time of Observable 1'


"""obs 2
mcmcs = [0.5497,0.5927,0.8715,1.8615,2.3618]
ecmcs = [0.5279,0.4971,0.4882,0.5102,0.4964]
mcmcerrors = [0.0189,0.0234,0.0449,0.1324,0.1838]
ecmcerrors = [0.0182,0.0100,0.0098,0.0144,0.0141]
observable = 'Integrated Autocorrelation time of Observable 2'
"""
"""
#Ac obs 3
mcmcs = [0.5046,0.4881,0.5009,0.5210,0.5060]
ecmcs = [0.4832,0.4922,0.4989,0.4891,0.4998]
observable = 'Integrated Autocorrelation time of Observable 3'
mcmcerrors = [0.0143,0.0098,0.0142,0.0180,0.0143]
ecmcerrors = [0.0097,0.0099,0.01,0.0098,0.0141]
"""

if len(mcmcs)==len(ecmcs):
	title = 'Plot of ' + observable + ' for 10000 sweeps (1000 Measurements)'
	plt.errorbar(xs[0:], mcmcs[0:], yerr=mcmcerrors[0:], color='red', label='Metropolis')
	plt.errorbar(xs[0:], ecmcs[0:], yerr=ecmcerrors[0:], color='blue', label='EventChain')
	plt.xlabel(r'Lattice size $L$')
	plt.ylabel(observable)
	plt.legend()
	#plt.ylim(0.4,0.6)
	plt.show()
	

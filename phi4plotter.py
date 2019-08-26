import matplotlib.pyplot as plt
import numpy as np
xs = [0.01,0.1,0.25,0.5,1.0]


#av action
#fill in these with data
mcmcs = [0.4997,
0.4997,
0.4997,
0.4999,
0.4998,
]
ecmcs = [0.4999,
0.5006,
0.5010,
0.5002,
0.5001,
]
expected = [0.5,0.5,0.5,0.5,0.5]
mcmcerrors = [0.0003,
0.0003,
0.0003,
0.0003,
0.0003,
]
ecmcerrors = [0.0004,
0.0004,
0.0004,
0.0004,
0.0005,
]
observable = 'Average Action per site'




"""

#obs1
mcmcs = [0.4993,
0.4963,
0.4843,
0.4529,
0.3728,
]
ecmcs =[0.4994,
0.4973,
0.4856,
0.4531,
0.3731,
]
expected = [0.4995,
0.4967,
0.4846,
0.4529,
0.3730,
]
mcmcerrors = [0.0004,
0.0003,
0.0003,
0.0003,
0.0002,
]
ecmcerrors = [0.0004,
0.0004,
0.0004,
0.0004,
0.0003,
]
observable = 'Observable 1'

"""


"""obs 2
mcmcs = [3.5578,
0.3400,
0.2467,
0.1883,
0.1270
]
ecmcs = [5.0521,
0.3311,
0.2471,
0.1885,
0.1269,
]
expected = [5.1828,
0.3321,
0.2467,
0.1884,
0.1270,
]
mcmcerrors = [1.7986,
0.0068,
0.0005,
0.0002,
0.0001,
]
ecmcerrors = [0.1876,
0.0011,
0.0004,
0.0002,
0.0001,
]
observable = 'Observable 2'
"""






"""obs 3
mcmcs = [3.2586,
0.0551,
0.0073,
0.0019,
0.0005,
]
ecmcs = [4.7523,
0.0483,
0.0078,
0.0019,
0.0005,
]
expected = [4.8778,
0.0485,
0.0076,
0.0018,
0.0004,
]
mcmcerrors = [1.7994,
0.0064,
0.0003,
0.0000,
0.0000,
]
ecmcerrors = [0.1876,
0.0009,
0.0002,
0.0000,
0.0000,
]
observable = 'Observable 3'

"""




"""obs 4
mcmcs = [51.7006,
43.2344,
20.2404,
6.9941,
1.9218,
]
ecmcs = [51.1054,
40.7643,
19.8255,
7.0338,
1.9079,
]
expected = [51.91043639,41.40043953,20.12544672,7.367486489,2.407481824]
mcmcerrors = [1.3106,
1.0624,
0.3335,
0.0794,
0.0195,
]
ecmcerrors = [0.5429,
0.4128,
0.1998,
0.0758,
0.0227,
]
observable = 'Observable 4'

"""



"""
#ac action
mcmcs = [0.4857,
0.5140,
0.4808,
0.5084,
0.4937,
]
ecmcs = [0.7339,
0.8070,
0.8616,
0.9413,
1.0897,
]
expected = []
mcmcerrors = [0.0139,
0.0251,
0.0193,
0.0203,
0.0140,
]
ecmcerrors = [0.0454,
0.0545,
0.0579,
0.0683,
0.0840,
]
observable = 'Integrated Autocorrelation time of Average Action'

"""




"""Ac obs 1
mcmcs = [0.6162,
0.5528,
0.5606,
0.5667,
0.5062,
]
ecmcs = [0.7303,
0.8069,
0.8570,
0.9449,
1.0838
]
expected = []
mcmcerrors = [0.0138,
0.0204,
0.0191,
0.0203,
0.0140,
]
ecmcerrors = [0.0452,
0.0545,
0.0576,
0.0685,
0.0836,
]
observable = 'Integrated Autocorrelation time of Observable 1'
"""

"""

Ac obs 2
mcmcs = [436.3843,
16.2528,
1.0695,
0.5271,
0.4869,
]
ecmcs = [1.9724,
0.5150,
0.6423,
0.7783,
1.0122,
]
expected = []
mcmcerrors = [188.2353,
3.6305,
0.0825,
0.0257,
0.0139,
]
ecmcerrors = [0.1896,
0.0205,
0.0357,
0.0479,
0.0730,
]
observable = 'Integrated Autocorrelation time of Observable 2'

"""




"""Ac obs 3
mcmcs = [436.6405,
18.3325,
1.7863,
0.6435,
0.4913,
]
ecmcs = [1.9715,
0.4969,
0.5104,
0.5852,
0.6022,
]
expected = []
mcmcerrors = [188.2419,
4.2854,
0.1655,
0.0358,
0.0197,
]
ecmcerrors = [0.1895,
0.0141,
0.0250,
0.0328,
0.0337,
]
observable = 'Integrated Autocorrelation time of Observable 3'

"""


"""Ac obs 4
mcmcs = [3.2278,
2.9272,
1.3635,
0.6372,
0.4956,
]
ecmcs = [0.5607,
0.5234,
0.5072,
0.5671,
0.5724,
]
expected = []
mcmcerrors = [0.3906,
0.3374,
0.1167,
0.0355,
0.0141,
]
ecmcerrors = [0.0272,
0.0256,
0.0248,
0.0318,
0.0418,
]
observable = 'Integrated Autocorrelation time of Observable 4'
"""


if len(mcmcs)==len(ecmcs):
	title = 'Plot of ' + observable + ' for 10000 sweeps (1000 Measurements)'
	plt.errorbar(xs[:], mcmcs[:], yerr=mcmcerrors[:], color='red', label='Metropolis')
	plt.errorbar(xs[:], ecmcs[:], yerr=ecmcerrors[:], color='blue', label='EventChain')
	if expected != []:
		plt.plot(xs[:], expected[:], color='black', label='Expected Numerically')
	plt.xlabel('Mass')
	plt.ylabel(observable)
	plt.ylim(0.49,0.51)
	plt.legend()
	#plt.ylim()
	plt.show()
	

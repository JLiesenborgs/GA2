import matplotlib.pyplot as plt
import sys
import pprint
import numpy as np

allValues = None
for l in open(sys.argv[1], "rt"):
    # is something like [ num ]
    l = l.split()[1:-1]
    if allValues is None:
        allValues = [ [] for i in range(len(l)) ]

    values = list(map(float, l))
    for i,x in enumerate(values):
        allValues[i].append(x)

muSigmaValues = [ (3,2), (5,1) ]
for dimValues,muSigma in zip(allValues, muSigmaValues):
    plt.hist(dimValues, bins=50, density=True)
    mu, sigma = muSigma
    x = np.linspace(-10,10,1000)
    px = 1/((2*np.pi)**0.5*sigma)*np.exp(-(x-mu)**2/(2*sigma**2))
    plt.plot(x, px)
    #plt.yscale('log')
    plt.show()


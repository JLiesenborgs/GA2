import numpy as np
import matplotlib.pyplot as plt

def pickIndexTournament(popSize, S):
    indices = np.random.randint(0,popSize-1,S)
    return indices.min()

def pickIndexRanked(popSize, b):
    x = np.random.random()
    y = 1.0 - x**(1.0/(1.0+b))
    idx = int(y*popSize)
    if idx >= popSize:
        idx = popSize-1
    return idx

def main():
    N = 128
    R = 1000000
    
    plt.figure(figsize=(12,6))
    idx = 0
    for pickFunction, fmtString, params in zip([pickIndexTournament, pickIndexRanked],
                                               ["Tournament size = {S}", "beta = {b}"],
                                               [ 
         #                                        [{"S": 1}, {"S": 2}, {"S": 3}, {"S": 4}],
         #                                        [{"b": 0}, {"b": 0.5}, {"b": 1.0}, {"b": 1.5}, {"b": 2.0}, {"b": 2.5}, {"b": 3}]
                                                [{"S": 4}],
                                                [{"b": 3}]
                                               ]):

        idx += 1
        #plt.subplot(1,2,idx)
        for p in params:

            counts = np.zeros(N, dtype=np.double)

            for _ in range(R):
                counts[pickFunction(N, **p)] += 1

            counts /= R
            print(counts)
            print(sum(counts))

            plt.plot(np.linspace(0.5, N-0.5, N), counts, '.', label=fmtString.format(**p))
            
        plt.legend()

    plt.show()


if __name__ == "__main__":
    main()

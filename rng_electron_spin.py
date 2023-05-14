from __future__ import division
import sys
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt

# this defines a function that is used to create an initial state of spin vectors
# this is used to start our monte carlo simulation of the ising model
def initialstate(N):   
    state = 2*np.random.randint(2, size=(N,N))-1
    return state

# this defines the monte carlo algoritm we will be using for this script
# it is a Markov chain Monte Carlo method, I am using severe help from a 
# website for information about this since I don't understand it entirely
# see sources for more information
def mcmove(config, beta):
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                s =  config[a, b]
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                cost = 2*s*nb
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                config[a, b] = s
    return config


# calculates the energy of the generated configuration
def calcEnergy(config):
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/4.


# calculates the magnetization of the configuration
def calcMag(config):
    mag = np.sum(config)
    return mag

# this class (also from online) helps define some ising model function
# see source information for more details
class Ising():
    ## monte carlo moves
    def mcmove(self, config, N, beta):
        for i in range(N):
            for j in range(N):            
                    a = np.random.randint(0, N)
                    b = np.random.randint(0, N)
                    s =  config[a, b]
                    nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                    cost = 2*s*nb
                    if cost < 0:	
                        s *= -1
                    elif rand() < np.exp(-cost*beta):
                        s *= -1
                    config[a, b] = s
        return config
    
    def simulate(self):   
        # initilizes the lattice simulation
        N, temp     = 64, .4
        config = 2*np.random.randint(2, size=(N,N))-1
        f = plt.figure(figsize=(15, 15), dpi=80);    
        self.configPlot(f, config, 0, N, 1);
        
        # want to create values at different times to be able to 
        # observe the evolution of the spins inside the lattice
        msrmnt = 1001
        for i in range(msrmnt):
            self.mcmove(config, N, 1.0/temp)
            if i == 1:       self.configPlot(f, config, i, N, 2);
            if i == 4:       self.configPlot(f, config, i, N, 3);
            if i == 32:      self.configPlot(f, config, i, N, 4);
            if i == 100:     self.configPlot(f, config, i, N, 5);
            if i == 1000:    self.configPlot(f, config, i, N, 6);
                    
    def configPlot(self, f, config, i, N, n_):
        # actually creates the plots mentioned before
        X, Y = np.meshgrid(range(N), range(N))
        sp =  f.add_subplot(3, 3, n_ )  
        plt.setp(sp.get_yticklabels(), visible=False)
        plt.setp(sp.get_xticklabels(), visible=False)      
        plt.pcolormesh(X, Y, config, cmap=plt.cm.RdBu);
        plt.title('Time=%d'%i); plt.axis('tight')    
    plt.show()

if __name__ == "__main__":
    if '-h' in sys.argv or '--help' in sys.argv:
        print("Usage: %s -seed [seed value] -Natoms [number of atoms in lattice] -Npoints [number of temperatures sampled]" % sys.argv[0])
        print
        sys.exit(1)

    ## default parameters
    # number of points of T used for calculation
    nt = 50

    # size of lattice, number of atoms (NxN)
    N = 10

    # steps for equilibration
    eqSteps = 1024
    mcSteps = 1024

    # default seed
    seed = 555      
    np.random.seed(seed)

    # these create a few arrays for calculation and plotting,
    # for temperature, energy, magnetization, heat capacity, and sucesptibility
    # (intensive values)
    T = np.linspace(1.53, 3.28, nt); 
    E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
    n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N) 
 
    # read the user-provided inputs from the command line (if there)
    if '-seed' in sys.argv:
        p = sys.argv.index('-seed')
        seed = int(sys.argv[p+1])
    if '-Natoms' in sys.argv:
        p = sys.argv.index('-Natoms')
        N = int(sys.argv[p+1])
    if '-Npoints' in sys.argv:
        p = sys.argv.index('-Npoints')
        nt = int(sys.argv[p+1])

    params = [nt, N, eqSteps, mcSteps, seed, n1, n2]
    
    x_vals = range(10000)
    y_vals = np.random.random(10000)

    
    # lets visualize our uniform distribution to check its uniformity
    plt.plot(x_vals, y_vals, '.')  
    plt.savefig('UniformDistribution.png')

    # this creates the equilibrium condition that the calculations will
    # need to satisfy
    for tt in range(nt):
        E1 = M1 = E2 = M2 = 0
        config = initialstate(N)
        iT=1.0/T[tt]; iT2=iT*iT;
    
        # does the monte carlo sim
        for i in range(eqSteps):
            mcmove(config, iT)

        # now we can calulcate the energy and magnetization
        for i in range(mcSteps):
            mcmove(config, iT)           
            Ene = calcEnergy(config)
            Mag = calcMag(config)   
            E1 = E1 + Ene
            M1 = M1 + Mag
            M2 = M2 + Mag*Mag 
            E2 = E2 + Ene*Ene

        E[tt] = n1*E1
        M[tt] = n1*M1
        C[tt] = (n1*E2 - n2*E1*E1)*iT2
        X[tt] = (n1*M2 - n2*M1*M1)*iT

    with open(r'energy.txt', 'w') as fp:
        for item in E:
            fp.write("%s\n" % item)

    with open(r'magnetization.txt', 'w') as fp:
        for item in M:
            fp.write("%s\n" % item)

    with open(r'heat_capacity.txt', 'w') as fp:
        for item in C:
            fp.write("%s\n" % item)

    with open(r'susceptibility.txt', 'w') as fp:
        for item in X:
            fp.write("%s\n" % item)

    with open(r'temperatures.txt', 'w') as fp:
        for item in T:
            fp.write("%s\n" % item)

    with open(r'parameters.txt', 'w') as fp:
        for item in params:
            fp.write("%s\n" % item)




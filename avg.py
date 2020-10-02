# Time-Averaging over trajectories run on different cores
import numpy as np
import matplotlib.pyplot as plt 

def plot_pops(tsteps, pops, dtc, ttraj, iwrite=1, savefig=0, showfig=1):
    if(iwrite): np.savetxt('pops.dat', np.vstack((tsteps, pops)).T, fmt='%18.10e')
    plt.xlabel('time (in a.u.)')
    plt.ylabel('ground state population')
    plt.ylim([0,1.1])
    plt.title('dtc= '+ str(dtc) +'; ntraj= ' + str(ttraj))
    plt.plot(tsteps,pops, marker='.')
    if (savefig): plt.savefig("pops.png")
    if (showfig): plt.show()

#===========================================================

p = np.loadtxt("input0.txt", dtype=int)
ncore = int(p[0])
ntraj = int(p[1])
ttraj = ncore*ntraj
dtc   = p[2]
ttime = p[3]
ntsteps = np.ceil(ttime/dtc).astype('int')

popsums = np.zeros(ntsteps)
for seed in range(1,ncore+1):
	fpath = './' + str(seed)+ '/popsum_'+str(ntraj)+'.dat'
	p = np.loadtxt(fpath, dtype='float')
	tsteps, psum = zip(*p)
	popsums = np.add(popsums, psum)
popsums_avg = popsums/ttraj

plot_pops(tsteps, popsums_avg, dtc, ttraj, iwrite=1, savefig=1, showfig=0)
#============================================================

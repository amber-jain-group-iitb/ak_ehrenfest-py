import numpy as np
import matplotlib.pyplot as plt

def ham_diab(q, imodel=11):

    if (imodel == 1):       
        # Tully Model- 1 : Simple avoided crossing 
        assert (len(q) == 1)
        nquant = 2

        # parameters
        A,B,C,D = 0.01, 1.6, 0.005, 1.0

        # Hamiltonian
        H = np.zeros([nquant, nquant])
        if (q >= 0):
            H[0,0] = A*(1 - np.exp(-B*q))
        if (q<0):
            H[0,0] = -A*(1 - np.exp(B*q))
        H[1,1] = - H[0,0]
        H[0,1] = C*np.exp(-D*(q**2))
        H[1,0] = H[0,1]

        # Hamiltonian derivative wrt q
        delq_H = np.zeros([nquant, nquant])
        if (q >= 0):
            delq_H[0,0] = A*B*np.exp(-B*q)
        if (q<0):
            delq_H[0,0] = A*B*np.exp(B*q)
        delq_H[1,1] = - delq_H[0,0]
        delq_H[0,1] = C*np.exp(-D*(q**2))
        delq_H[1,0] = delq_H[0,1]   

        # vibrionic coupling : non-adiabatic (derivative) coupling vector
        dc = np.zeros(H.shape, dtype=float)
        z = (H[0,0] - H[1,1])/(H[0,1]*2)
        dc[1,2] = 1/(2*(1+z**2))*( ( 1/(2*H[0,1]**2)) * 
                    ( H[0,1]*(delq_H[0,0]-delq_H[1,1]) - delq_H[0,1]*(H[0,0]-H[1,1]) ) ) 
        dc[2,1]= -dc[1,2]

    if (imodel == 11):      
        # Spin boson : 2 level system coupled to a harmonic bath
        # assert (len(q) == 1)
        nquant = 2

        # parameters
        # en    = 120.0
        # vcoup = 87.7
        # omega = 1060.0
        # g     = 2.61E-3
        # mn    = 1836.0

        # Hamiltonian
        H  = np.zeros([nquant, nquant], dtype=float)
        H =H+ np.array([[-en/2, vcoup],
                       [vcoup, en/2]])
        # harmonic potential
        H =H+ (0.5* mass * omega**2 * q**2)*np.eye(nquant, dtype=float)
        # nuclear kinetic energy
        # H =H+ 0.5* mass * q_dot**2
        # coupling between two level system and bath
        H =H+ np.array([[g*q, 0.0],
                        [0.0, -g*q]])

        # Hamiltonian derivative
        delq_H  = np.zeros([nquant, nquant])
        delq_H =delq_H+ (mass * omega**2 * q)*np.eye(nquant, dtype=float)
        # delq_H =delq_H+ force 
        delq_H =delq_H+ np.array([[g, 0], 
                                 [0, -g]])      
    
    # dictionary : ham_model['H', 'dH', 'dc']   # [nquant x nquant]*npart
    return H, delq_H#, dc

def d_ij(phi, E, delq_H):
    # vibrionic coupling : non-adiabatic (derivative) coupling
    dc = np.zeros(H.shape, dtype=float)
    for i in range(nquant):
        for j in range(i+1, nquant):
            dc[i,j] = (phi[i].T*delq_H*phi[j])/(E[i]-E[j])
        for j in range(0,i):
            dc[i,j] = -dc[j,i]
    return dc

def vdotd(q_dot):
    pass

def ham_adiab(H, q):
    eval, evec = np.linalg.eigh(ham_diab(q, imodel=11))
    return eval, evec

# Quantum Dynamics
def tdse_cdot(ci, H):
    iota = complex(0,-1)
    # return iota/hbar*np.dot(H,ci) - vdotd()
    # return  -vdotd()                      # adiabatic basis
    return iota/hbar*np.dot(H,ci)       # dibatic basis

def evolve_quantum(ci, H_tn_1, H_tn_4, dtc):
    # Runge-Kutta 4th order integrator
    k1 = tdse_cdot(ci,H_tn_1)
    
    H_tn_23 = (H_tn_1 + H_tn_4)/2
    k2 = tdse_cdot(ci+dtc*k1/2, H_tn_23) 
    k3 = tdse_cdot(ci+dtc*k2/2, H_tn_23)  
    
    k4 = tdse_cdot(ci+dtc*k3, H_tn_4)

    ci = ci + dtc*(k1+2*k2+2*k3+k4)/6
    return ci

# Gaussian sampling over trajectory variables
def gauss_rand_n(seed=None, mu=0, sig=1):
    np.random.seed(seed)
    r1,r2   = np.random.rand(2)
    r       = np.sqrt(-2*np.log(r1))*np.cos(2*np.pi*r2)
    rnd     = mu + sig*r
    return rnd

# Classical Dynamics
def compute_force(ci, delH_diab):
    force = - np.dot(np.dot(ci.conj().T, delH_diab), ci)
    return np.real(force)
    
def stochasticforce(dtc, imodel = 11):
    # stochastic force : Langevin  dynamics??
    assert(ifriction == 1)
    if (imodel == 11): 
        # omega = 1060.0 #
        # eta   = 10*omega

        # kb    = 1.3806503E-23
        # tempK = 77.0
        # beta  = 1/(kb*tempK)
        # mass

        gamma_B = eta
        gdt = gamma_B*dtc  #

        npart = len(mass)
        delr, delv = np.zeros(npart), np.zeros(npart)
        for i in range(npart):    # npart    
            sig_r=dtc*np.sqrt(1/(beta*mass[i]) *1/gdt*(2-1/gdt*(3-4*np.exp(-gdt)+np.exp(-2*gdt))))
            sig_v=np.sqrt(1/(beta*mass[i])*(1-np.exp(-2*gdt)))
            sig_rv=(dtc/(beta*mass[i])* 1/gdt *(1-np.exp(-gdt))**2)/(sig_r*sig_v)  # correlation coeffecient
            
            # call gaussian_random_number(rnd1)
            # call gaussian_random_number(rnd2)
            rnd1, rnd2  = gauss_rand_n(), gauss_rand_n()
            delr[i]     = sig_r*rnd1
            delv[i]     = sig_v*(sig_rv*rnd1+np.sqrt(1-sig_rv**2)*rnd2)
        #  delr=delr-sum(delr)/dfloat(nclass)
        #  delv=delv-sum(delv)/dfloat(nclass     
    return delr, delv

def evolve_classical(q, q_dot, force, dtc, ifriction):
    # Velocity Verlet
    if (ifriction == 0):
        #force = -delq_E  # ham_model(q)['H','dH']
        acc = force/mass
        q = q + q_dot*dtc + 0.5*acc*(dtc**2)

        H, delq_H = ham_diab(q)
        force = compute_force(ci, delq_H)

        acc = 0.5*(acc + force/mass)
        q_dot = q_dot + acc*dtc

    elif (ifriction == 1):
        gamma_B = eta
        gdt = gamma_B*dtc
        c0 =np.exp(-gdt)
        c1 =1/gdt*(1-c0)
        c2 =1/gdt*(1-c1)

        acc = force/mass
        delr, delv = stochasticforce(dtc)
        q = q+c1*dtc*q_dot+ c2*dtc**2*acc+ delr

        H, delq_H = ham_diab(q)
        force = compute_force(ci, delq_H)

        acc = (c1-c2)*acc + force/mass
        q_dot = c0*q_dot+ acc*dtc +delv

    return q, q_dot, H, force

# Population readout
def population(ci, istate):
    return np.abs(ci[istate])**2

def check_ci(ci):
    return np.sum(np.abs(ci)**2)

def check_energy(ci, H_diab):
    #assert(ifriction == 0)
    #?? return 0.5*mass*q_dot**2 + 0.5*mass*omega**2*q**2
    return np.real( np.dot(np.dot(ci.conj().T, H_diab), ci) ) # +0.5*mass*q_dot**2

#..........................................................#
# plotting functions
#..........................................................#

def plot_surf(qs):
    fig = plt.figure(1)
    H00, H11, H01 = [], [], []
    for q in qs:
        H, delq_H = ham_diab(q)
        Hplot = H  # H, delq_H
        H00.append(Hplot[0,0])
        H11.append(Hplot[1,1])
        H01.append(Hplot[0,1])
    plt.plot(qs, H00, color='b', marker='.', label='00')
    plt.plot(qs, H11, color='g', marker='.', label='11')
    plt.plot(qs, H01, color='r', marker='.', label='01')
    plt.xlabel("q (in a.u.)")
    plt.ylabel("energy")
    plt.legend()
    plt.show()

def plot_pops(tsteps, pops, dtc, ntraj, savefig=0, showfig=1):
    fig = plt.figure(2)
    plt.xlabel('time (in a.u.)')
    plt.ylabel('ground state population')
    plt.ylim([0,1.1])
    plt.title('dtc= '+ str(dtc) +'; ntraj= ' + str(ntraj))
    plt.plot(tsteps,pops, marker='.')
    if (savefig): plt.savefig("pops.png")
    if (showfig): plt.show()

def plot_tot_energy(tsteps, tot_energy, itraj=1, savefig=0, showfig=0):
    fig = plt.figure(3)
    plt.xlabel('time (in a.u.)')
    plt.ylabel("energy")
    #plt.ylim([-0.1,1])
    plt.plot(tsteps, tot_energy[itraj], marker='.')
    if (savefig): plt.savefig("eng.png")
    if (showfig): plt.show()

def plot_sum_ci2(tsteps, sum_ci_sq, itraj=1, savefig=0, showfig=0):
    fig = plt.figure(4)
    plt.xlabel('time (in a.u.)')
    plt.ylabel("sum_ci")
    #plt.ylim([0,2])
    plt.plot(tsteps,sum_ci_sq[itraj], marker='.')
    if (savefig): plt.savefig("sum_ci.png")
    if (showfig): plt.show()
#..........................................................#
#..........................................................#


# unit conversions
hbar   = 1.0
clight = 2.99792458E10          # cm/sec
sec2au = 4.13414E16             # sec2au 
omg2au = clight*2*np.pi/sec2au	# cm_1 to au
eng2au = 4.55634E-6             # cm_1 to au
tk2au  = 3.16681E-6				# K to au
fs2au  = 40.0                   # 41.34137au
# parameters
en    = 120.0*eng2au
vcoup = 87.7*eng2au
kb    = 1.0
omega = 1060.0*omg2au
eta   = 10*omega
g     = 2.61E-3 # a.u.
tempK = 77.0*tk2au
beta  = 1/(kb*tempK)

npart  = 1
nquant = 2
mass  = 1836.0*np.ones(npart)

sig_q     = 1/np.sqrt(beta*mass*omega**2)
sig_q_dot = 1/np.sqrt(beta*mass)

# Read parameters
param_ll = np.loadtxt('input0.txt', dtype=int)
iflow = param_ll[0]
# if (iflow == 1):
# 	param_ll = np.loadtxt('input.txt', dtype=int)

#----------------------------------------------------------#
# PES
#----------------------------------------------------------#
qs, qstep = np.linspace(-0.25, 0.25, 1001, retstep=True, endpoint=True)
if (param_ll[2] == 0):  	# ntraj
    print(qstep, qs)
    plot_surf(qs)
#----------------------------------------------------------#

#----------------------------------------------------------#
# Dynamics
#----------------------------------------------------------#

iseed       = param_ll[-1]
ifriction   = param_ll[5]

ntraj       = int(param_ll[2])
dtc         = param_ll[3]   #param_ll[2]
ttime       = param_ll[4]   #param_ll[3]*fs2au
nsteps = np.ceil(ttime/dtc).astype('int')

istate	= param_ll[6]
iprint	= param_ll[7]
istr_sci= param_ll[8]
istr_en	= param_ll[9]

sum_ci_sq   = np.zeros([ntraj, nsteps], dtype=float)
tot_energy  = np.zeros([ntraj, nsteps], dtype=float)
pop         = np.zeros([ntraj, nsteps], dtype=float)
for nt in range(ntraj):
    rnd1, rnd2  = gauss_rand_n(), gauss_rand_n()
    q           = (-g/(mass*omega**2)+rnd1*sig_q)   #*np.ones(npart, dtype=float)
    q_dot       = (rnd2*sig_q_dot)*np.ones(npart)   #*np.ones(npart, dtype=float)
    H, delq_H   = ham_diab(q)                       # [nquant x nquant]*npart
    
    ci          = np.array([complex(1,0), complex(0,0)]) # nquant*npart
    force       = compute_force(ci, delq_H)         #*npart
    for ns in range(0,nsteps):
        if (iprint==1): print(nt, ns, check_ci(ci), check_energy(ci, H))	# print(q,q_dot)
        if (istr_sci): sum_ci_sq[nt,ns]  = check_ci(ci)
        if (istr_en) : tot_energy[nt,ns] = check_energy(ci, H)
        pop[nt,ns]   = population(ci, istate)
        # evolve trajectory
        q, q_dot, H_new, force  = evolve_classical(q, q_dot, force, dtc, ifriction)
        ci                      = evolve_quantum(ci, H, H_new, dtc)
        H                       = H_new
        #if np.abs(check_ci(ci)-1.0) > 0.001 : break
    if (iprint==2): print(nt, ns, check_ci(ci), check_energy(ci, H))
    popsum  = np.sum(pop,axis=0)
    pops    = popsum/ntraj
#----------------------------------------------------------#
tsteps= [dtc*i for i in range(0,nsteps)]
if (iflow == 0):
	# average population
	np.savetxt('popsum_'+str(ntraj)+'.dat', np.vstack((tsteps, pops)).T, fmt='%18.10e')
if (iflow == 1):
	# population_sum
	np.savetxt('popsum_'+str(ntraj)+'.dat', np.vstack((tsteps, popsum)).T, fmt='%18.10e')

# plots: population; total_energy, |ci|^2 for 2nd trajectory
plot_pops(tsteps,pops, dtc, ntraj, savefig=1, showfig=0)
if (istr_sci): plot_sum_ci2(tsteps, sum_ci_sq, itraj=1, savefig=1, showfig=0)
if (istr_en) : plot_tot_energy(tsteps, tot_energy, itraj=1, savefig=1, showfig=0)
#----------------------------------------------------------#



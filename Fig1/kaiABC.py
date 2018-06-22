"""
System has N "KaiC hexamers", these molecules go through a cycle from
C_bot -> C1_up -> C2_up -> ... C(m-1)_up -> C_top -> C(m-1)_down -> C(m-2)_down ->
C_bot
 
m parameterizes the number of steps in the cycle

On the down part of the cycle, transitions are only "downhill" towards
C_bot, and proceed with a rate k* = 2mk.
 
On the up part of the cycle, the activator A ("KaiA") can bind to give
A*C_bot, A*C1_up, etc.  binding occurs with 2nd order rate constant
k_on_up, and dissociation with k_off_up.  Assume these are the same for
the entire up half of the cycle.
 
A*C complexes go "up" with rate k*.  C_up complexes not bound to A fall
back "down (e.g. C2_up -> C1_up) with rate k*.
 
On the down part of the cycle, A can bind to give C_top*A, C(m-1)_down*A,
etc. with 2nd order rate constant k_on_down.  For now, assume A does not
dissociate until the cycle is completed (at C_bot).
 
Thus, m parameterizes the model topology. N specifies the scale of the
simulation.
 
k,  k_on_up, k_off_up, k_on_down, are kinetic constants.

Justin Chew
May 23, 2017
Rust Lab
*modified by Weerapat Pittayakanchit
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.interpolate as interpolate
import sys

################################
## Function definitions below ##
################################

def generate_species(pickle_file=None):
    """
    Generates all possible molecular species in the simulation, where the type of
    molecule is specified as a tuple (p, a, n).

    p refers to whether the KaiC hexamer is on the up cycle (phosphorylating, p = 1)
    or down cycle (dephosphorylating, p = 0).

    a refers to whether KaiC hexamer is bound to KaiA (a = 1 for bound, a = 0 for unbound),
    and whether it is sequestering or activating depends on whether p = 0 or p = 1, respectively.
    a = 2 represents a form of KaiC, C*, that can phosphorylate in the absence of KaiA (i.e. the
    A loop tails are exposed).

    n refers to the number of phosphorylations on the hexamer and can range from 0 to m.

    The paramter pickle_file species what pickle file the function should take in
    to be used to initialize the amounts of each Kai species to some specific amount.
    """

    KaiC_species = {}

    if not pickle_file:
        for p in range(2):
            for a in range(3):
                for n in range(m):
                    if p:
                        KaiC_species[(p, a, n)] = 0
                    else:
                        if a < 2:
                            KaiC_species[(p, a, n+1)] = 0
        KaiC_species[(1, 0, 0)] = KaiC_init

    return KaiC_species


def generate_rxns(species):
    """
    Given a particular molecular species, this will return all the possible reactions
    that destroy this molecule.  Only counting reactions that destroy the molecule avoids
    double counting reactions, which works because there are no reactions that create or
    destroy molecules.

    Reactions are structured as a tuple (species1, species2, type), where species1 is the
    reactant, species2 is the product, and type refers to the reaction type:
    type = 0 -> phosphorylation/dephosphorylation
           1 -> KaiA binding on up cycle
           2 -> KaiA unbinding on up cycle
           3 -> KaiA sequestering on down cycle
           4 -> phosphorylating to Cm and unbinding KaiA
           5 -> dephosphorylating to C0 and unbinding KaiA
           6 -> transition to C* on up cycle
           7 -> transition from C* on up cycle
    """

    rxn_list = []

    if species[0]:   # add rxns for up cycle KaiC
        if not species[1]:  # add rxns for unbound KaiC
            if species[2]:
                rxn_list.append((species, (1, 0, species[2]-1), 0))   # dephosphorylate
            rxn_list.append((species, (1, 1, species[2]), 1)) # bind to KaiA
            rxn_list.append((species, (1, 2, species[2]), 6))   # transition to C*
        if species[1]:   # add rxns for KaiC bound to KaiA and for C*
            if species[2] < m-1:
                rxn_list.append((species, (1, species[1], species[2]+1), 0))   # phosphorylate
            if species[2] == m-1:
                rxn_list.append((species, (0, 0, m), 4))  # phosphorylate to m phos state and unbind KaiA
            rxn_list.append((species, (1, 0, species[2]), 2)) # unbind KaiA or transition back to C

    if not species[0]:   # add rxns for down cycle KaiC
        if species[2] == 1: # add rxns specific for m = 1 down cycle KaiC bound/unbound
            if species[1]:
                rxn_list.append((species, (1, 0, 0), 5)) # dephos to C0 and unbind KaiA
            if not species[1]:
                rxn_list.append((species, (1, 0, 0), 0))  # dephos to C0
        if species[2] > 1: # add rxns specific for all m > 1 down cycle KaiC bound/unbound
            rxn_list.append((species, (0, species[1], species[2]-1), 0))  # dephos by one
        if not species[1]:  # add rxns specific for unbound KaiC
            if species[2] == m: # add rxn for totally phos KaiC back to m-1 KaiC on up cycle
                rxn_list.append((species, (1, 0, m-1), 0))
            rxn_list.append((species, (0, 1, species[2]), 3)) # sequester KaiA

    return rxn_list


def calc_rates():
    """
    Returns the reaction rate for all reactions based on the
    type of reaction, specified in the third element of the rxn tuple:

    type = 0 -> phosphorylation/dephosphorylation
           1 -> KaiA binding on up cycle
           2 -> KaiA unbinding on up cycle
           3 -> KaiA sequestering on down cycle
           4 -> phosphorylating to Cm and unbinding KaiA
           5 -> dephosphorylating to C0 and unbinding KaiA
           6 -> transition to C* on up cycle
           7 -> transition from C* on up cycle
    """

    rxn_rates = []
    
    for rxn in rxns:
        rxn_type = rxn[2]
        if rxn_type == 0 or rxn_type == 4 or rxn_type == 5:
            rate = kstar
        if rxn_type == 1:
            if A_init:
                rate = k_on_up*conc*(A/A_init)*eta*ATP
            else:
                rate = 0.0
        if rxn_type == 2 or rxn_type == 7:
            rate = k_off_up
        if rxn_type == 3:
            if A_init:
                rate = k_on_down*conc*(A/A_init)
            else:
                rate = 0.0
        if rxn_type == 6:
            rate = k_on_up*(1-eta)*ATP

        rxn_rates.append(rate*KaiC_species[rxn[0]])

    return rxn_rates


def calc_P_KaiC(KaiC_species):
    """
    Returns the amount of phosphorylated KaiC.
    """

    return sum([species[2]*KaiC_species[species] for species in KaiC_species])


def calc_KaiBC_complexes(KaiC_species):
    """
    Returns the amount of KaiC hexamers that are sequestering KaiA.
    """

    return sum([KaiC_species[species] for species in KaiC_species if not species[0] and species[1]])


def record_data():
    """
    Records current datapoint into the global data list for output.
    """

    global data

    temp_data = []
    temp_data.append(t)
    temp_data.append(calc_P_KaiC(KaiC_species))
    temp_data.append(calc_KaiBC_complexes(KaiC_species))
    # nonseq KaiA temp_data.append(A_init - stoich*sum([species[1]*KaiC_species[species] for species in KaiC_species if species[0] == 0]))
    temp_data.append(A_init - stoich*sum([species[1]*KaiC_species[species] for species in KaiC_species if species[0] == 0]) - calc_free_KaiA())

    if sum([KaiC_species[x] for x in KaiC_species if x[0] == 1]):
        avg_phos_KaiC = 1.0/sum([KaiC_species[x] for x in KaiC_species if x[0] == 1])*sum([x[2]*KaiC_species[x] for x in KaiC_species if x[0] == 1])
    else:
        avg_phos_KaiC = np.nan
    temp_data.append(avg_phos_KaiC)

    if sum([KaiC_species[x] for x in KaiC_species if x[0] == 0]):
        avg_phos_KaiC_d = 1.0/sum([KaiC_species[x] for x in KaiC_species if x[0] == 0])*sum([x[2]*KaiC_species[x] for x in KaiC_species if x[0] == 0])
    else:
        avg_phos_KaiC_d = np.nan
    temp_data.append(avg_phos_KaiC_d)
    temp_data.append(sum([KaiC_species[x] for x in KaiC_species if x[0]==0 and x[1]==0]))
    temp_data.append(ATP)
    temp_data.append(1.0*sum([KaiC_species[x] for x in KaiC_species if x[0]==1])/sum(KaiC_species.values()))
    temp_data.append(1.0*sum([KaiC_species[x] for x in KaiC_species if x[0]==0])/sum(KaiC_species.values()))
    temp_data.append(sum([x[2]*KaiC_species[x] for x in KaiC_species if x[0]==1])/(m*KaiC_init/2))
    temp_data.append(sum([x[2]*KaiC_species[x] for x in KaiC_species if x[0]==0])/(m*KaiC_init/2))

    data.append(temp_data)


def calc_free_KaiA():
    """
    Returns the amount of free KaiA.
    """

    bound_A = sum([species[1]*KaiC_species[species] for species in KaiC_species if species[0] == 1 and species[1] == 1])
    sequest_A = stoich*sum([species[1]*KaiC_species[species] for species in KaiC_species if species[0] == 0])
    A = max(0, A_init - (bound_A + sequest_A))

    return A


def RKfunction(rates, KaiC_species):
    """
    Calculates each iteration of Runga Kutta with the functions defined
    herein (e.g. the systems of ODEs that drive the system). Takes as
    arguments a list of rates for each reaction as well as the
    dictionary of all possible KaiC species and returns a dictionary of
    dspecies/dt for each species.
    """
    
    # output contains the ODE functions iterated once on the input variables
    output = {}
    
    for species in KaiC_species:
        
        dspecies = 0.0
        
        for reaction in reaction_dict_f[species]:
            dspecies += rates[reaction[1]]
        for reaction in reaction_dict_r[species]:
            dspecies -= rates[reaction[1]]
        
        output[species] = dspecies
    
    return output


def compile_rxn_dicts():
    """
    Compiles the forward and reverse reaction dictionaries ahead of time
    for the RK4 algorithm in the deterministic version of the simulation.

    The forward and reverse reaction dictionaries hold the reactions that
    create and destroy a particular species, respectively.

    Each entry in the dictionary also holds the respective index of the rxn
    in the rxns list so that the appropriate rate can be pulled from the
    rates list.
    """

    reaction_dict_f = {}
    reaction_dict_r = {}

    for species in KaiC_species:
        reaction_dict_f[species] = []
        reaction_dict_r[species] = []

    for index, rxn in enumerate(rxns):
        reactant = rxn[0]
        product = rxn[1]
        reaction_dict_f[product].append((rxn, index))
        reaction_dict_r[reactant].append((rxn, index))

    return reaction_dict_f, reaction_dict_r


def get_ATP_level(t):
    """
    Returns ATP level based on current time.
    """

    # if env_noise and LD_flag:
    #     return noise_func(t)

    # else:
        #if t > LD_start and t < (LD_start + LD_period*LD_cycles) and LD_flag:
    OSCILLATION_Flag = True
    if OSCILLATION_Flag:
        if ((t-LD_start) % LD_period) < (LD_period*day_frac):
            if env_noise==1 and LD_flag:
                return noise_func(t)
            else:
                return ATP_range[1]
        else:
            return ATP_range[0]
    else:
        return ATP_range[1]


def generate_env_noise():
    """
    Generates a timecourse of environmental input noise where the timing of
    light and dark are deterministic but the actual light levels are stochastic.
    The light levels are drawn from a normal distribution centered on the high or
    low ATP level specified in the parameters and where the transition times are
    drawn from an exponential distribution with a given waiting time.

    The function returns an interpolated function that returns an ATP level given
    the time (using zero order interpolation).
    """

    k_change = 1/0.75    # rate is 1/h
    light_std = 0.25

    t = 0.0
    end_t = endtime

    t_events = []
    light_events = []

    DL_transitions = np.linspace(LD_start, LD_start+LD_cycles*LD_period, LD_cycles+1)
    LD_transitions = np.linspace(LD_start+LD_period*day_frac, LD_start+LD_period*day_frac+(LD_cycles-1)*LD_period, LD_cycles)

    # Keep track of which transition to test for
    LD_index = 0
    DL_index = 0

    ATP = ATP_range[1]

    # Flag used to end one more point past the end time for interpolator to cover whole time range
    run_once_more = True

    while t < end_t or run_once_more:

        # Clip light level to between 0 and 1
        # light_level = np.clip(np.random.normal(ATP, light_std), ATP_range[0], ATP_range[1])
        light_level = np.random.uniform(ATP_range[0], ATP_range[1])

        t_events.append(t)
        light_events.append(light_level)

        # print(str(correlation_dt/dt))
        noise_dt = np.random.exponential(1/correlation_dt)

        if LD_index <= LD_cycles - 1:
            if t + noise_dt > LD_transitions[LD_index]:
                noise_dt = LD_transitions[LD_index] - t
                ATP = ATP_range[0]
                LD_index += 1
        if DL_index <= LD_cycles:
            if t + noise_dt > DL_transitions[DL_index]:
                noise_dt = DL_transitions[DL_index] - t
                ATP = ATP_range[1]
                DL_index += 1

        if not t < end_t:
            run_once_more = False

        t += noise_dt

    interp_func = interpolate.interp1d(t_events, light_events, kind="zero")

    return interp_func

def ladder_interpolate(t, y, t_interpolate):
    # Assume that t_interpolate[0] > t[0] and t_interpolate[end] < t[end]
    y_interpolate = np.zeros(len(t_interpolate))
    index = 0
    for i in range(len(t_interpolate)):
        while t[index] < t_interpolate[i]:
            index = index + 1

        if index == 0:
            y_interpolate[i] = y[0]
        else:
            y_interpolate[i] = y[index - 1]

    return y_interpolate

def find_mi(time, x, y):
    MI_dt = 0.5
    t_bins = np.arange(0, 24, MI_dt)

    dx = 0.01
    dy = 0.01
    x_bins = np.arange(min(x), max(x) + dx, dx)
    y_bins = np.arange(min(y), max(y) + dx, dy)
    # print('min = ' + str(min(x)) + ', max = ', str(max(x)))
    # print(x_bins)
    # print(y_bins)

    count_xy, _, _ = np.histogram2d(x, y, bins=(x_bins, y_bins))
    Pxy = count_xy[np.nonzero(count_xy)]/len(x)

    Hxy = - sum(Pxy * np.log2(Pxy))

    temp_x = [[] for _ in range(len(t_bins))]
    temp_y = [[] for _ in range(len(t_bins))]
    HxyGivenT = 0

    for i in range(len(time)):
        time_of_day = np.mod(time[i], 24)
        MI_t_index  = int(np.floor(time_of_day/MI_dt))
        temp_x[MI_t_index].append(x[i])
        temp_y[MI_t_index].append(y[i])

    for i in range(len(t_bins)):
        count_xyt, _, _ = np.histogram2d(temp_x[i], temp_y[i], bins=(x_bins, y_bins))
        PxyGivenT = count_xyt[np.nonzero(count_xyt)]/len(temp_x[i])
        # print np.nonzero(count_xyt)
        # print PxyGivenT
        # print(sum(PxyGivenT))
        # print 1.0*len(temp_x[i])/len(time)
        HxyGivenT += - (1.0*len(temp_x[i])/len(time)) * sum(PxyGivenT * np.log2(PxyGivenT))
        # print('HxyGivenT = ' + str(HxyGivenT))

    MI = Hxy - HxyGivenT
    return MI

####################################
## Set simulation parameters here ##
####################################

if len(sys.argv) < 7:
    print('Program needs 7 parameters:')
    print('python kaiABC.py [ATP_min] [ATP_max] [eta] [sfactor] [stoch_flag] [env_noise]')
    print('Try:python kaiABC.py 0.2 0.8 1 1 0 0')
    print('#copy number = sfactor*1200')
    exit(1)

ATP1       = float(sys.argv[1])
ATP2       = float(sys.argv[2])
eta        = float(sys.argv[3])# contribution of KaiA to phosphorylation vs. just ATP
sfactor    = float(sys.argv[4])# copy number = sfactor*1200
stoch_flag = int(sys.argv[5])  # 0 is determinstic, 1 is stochastic
env_noise  = int(sys.argv[6])  # 0 means clean square wave, 1 means noisy square wave

if stoch_flag:
    print('stoch_flag = true')
if env_noise:
    print('env_noise = true')
sfactor = 1   # scale factor to adjust total number of molecules in system
m = 18          # number of steps to reach full phosphorylation

k = 0.2           # sets the characteristic phos/dephos rates (adjusted later for m)
k_on_up = 1.0     # rate at which KaiC binds KaiA during phosphorylation half of cycle
k_off_up = 0.1  # rate at which KaiC-KaiA disassociates during phos half of cycle
k_on_down = 5e2 # rate at which KaiC binds KaiA during dephosphorylation

conc = 2.0      # concentration of KaiA in uM (for both deterministic and stochastic)
stoich = 6      # the number of KaiA dimers inhibited per hexamer
AC_ratio = 1.0  # the relative stoichiometry of KaiA dimers to KaiC hexamers

LD_flag = True
LD_cycles = 4
LD_start = 0
LD_period = 24
day_frac = 0.5
ATP_range = (ATP1, ATP2)

starttime = 24*10.0
endtime = 24*20.0  # the total amount of time to simulate
correlation_dt = 1
dt = 0.01     # timestep for deterministic simulation

constant_vol = False

#####################################
## Initialize simulation variables ##
#####################################

if constant_vol:
    conc = conc * sfactor

# Scale reaction rates to give 24 hr oscillations
# adjust = 1.0/4.08
adjust = 1.0/4.055
k *= adjust
k_on_up *= adjust
k_off_up *= adjust
k_on_down *= adjust

# Initialize amounts of KaiC and KaiA
if stoch_flag == 0:
    KaiC_init = conc / AC_ratio
    A_init = conc
elif stoch_flag == 1:
    KaiC_init = round(1200 * sfactor)
    A_init = round(KaiC_init)

kstar = k*2*m       # scaled by 2*m, or the total number of steps in the simulation

# Set up dictionary entries for all species of KaiC hexamer
KaiC_species = generate_species()

# Set up all possible reactions
rxns = []
for species in KaiC_species:
    rxns.extend(generate_rxns(species))

# Compile all forward and reverse reactions for all species
reaction_dict_f, reaction_dict_r = compile_rxn_dicts()

# Initialize reaction rate list (order corresponds to rxns list)
rates = []

################################
## Execute and run simulation ##
################################

# Initalize noise in environment
if env_noise == 1:
    noise_func = generate_env_noise()

t = 0
ATP = ATP_range[1]

# Initialize data variable (eventually 2d np array) to hold sim data
data = []
record_data()

startTime = datetime.now()

if stoch_flag == 0:

    while t < endtime:

        ATP = get_ATP_level(t)

        A = calc_free_KaiA()
        rates = calc_rates()

        # Initialize list to hold interative versions of the KaiC_species
        # dictionary where indices 1 through 4 hold the four Runge Kutta
        # calculation states in the algorithm
        d = []
        for i in range(6):
            d.append(None)    # dicts hold all of the concentrations for each species

        # initialize first data point concentrations
        d[0] = {species: KaiC_species[species] for species in KaiC_species}

        # begin RK calculations http://www.myphysicslab.com/runge_kutta.html
        d[1] = RKfunction(rates, d[0])
        d[2] = RKfunction(rates, {key: d[0][key] + dt/2 * d[1][key] for key in d[1].keys()})
        d[3] = RKfunction(rates, {key: d[0][key] + dt/2 * d[2][key] for key in d[2].keys()})
        d[4] = RKfunction(rates, {key: d[0][key] + dt * d[3][key] for key in d[3].keys()})

        # calculate final differences in concentrations
        d[5] = {key: dt/6 * (d[1][key] + 2 * d[2][key] + 2 * d[3][key] + d[4][key]) for key in d[0].keys()}

        # update all values with RK values
        for species in KaiC_species:
            KaiC_species[species] += d[5][species]
            
        t += dt

        record_data()

if stoch_flag == 1:

    while t < endtime:
        
        ATP = get_ATP_level(t)

        A = calc_free_KaiA()
        rates = calc_rates()

        r_tot = sum(rates)
        simulate_dt = np.random.exponential(1/r_tot)

        selector = 0.0
        pick = np.random.random()*r_tot
        for index, rxn in enumerate(rxns):
            if pick < (selector + rates[index]):
                KaiC_species[rxn[0]] -= 1
                KaiC_species[rxn[1]] += 1
                break
            else:
                selector += rates[index]

        t += simulate_dt

        record_data()

###################################
## Plot and write data to output ##
###################################

endTime = datetime.now()

data = np.array(data)

header = ",".join(["time", "PKaiC-up", "PKaiC-dn"])
base_name = "test_output.csv"
special_name = str(eta)
initial_index_of_data = [n for n, i in enumerate(data[:,0]) if i > starttime][0]
#last_index_of_data    = [n for n, i in enumerate(data[:,0]) if i > endtime - 24][0]
np.savetxt("eta =  " +str(eta) + " ATP = " + str(ATP_range) + " copy_number = " + str(sfactor *1200) + ", stochastic = " + str(stoch_flag) + ", env_noise = " + str(env_noise) + ".csv", data[initial_index_of_data:,[0,10,11]], delimiter=",", header=header)

fig = plt.figure(figsize=(8,10))
ax = plt.subplot(111)
ax.plot(data[:,0], data[:,7], label="ATP", color="gray")
ax.plot(data[:,0], data[:,11], label="PKaiC dn", color="green")
ax.set_xlabel("Time (h)")

# ax = plt.subplot(711)
# ax.plot(data[:,0], data[:,1], label="P-KaiC", color="black")
# ax.legend()
# ax2 = plt.subplot(712)
# ax2.plot(data[:,0], data[:,2], label="KaiBC", color="gray")
# ax2.legend()
# ax3 = plt.subplot(713)
# ax3.plot(data[:,0], data[:,3], label="Phos. KaiA", color="red")
# ax3.legend()
# ax4 = plt.subplot(714)
# ax4.plot(data[:,0], data[:,4], label="Avg up cycle m", color="blue")
# ax4.plot(data[:,0], data[:,5], label="Avg dn cycle m", color="green")
# ax4.set_ylim(0, m)
# ax4.legend()
# ax5 = plt.subplot(715)
# ax5.plot(data[:,0], data[:,6], label="Unbd dn KaiC", color="black")
# ax5.legend()
# ax6 = plt.subplot(716)
# ax6.plot(data[:,0], data[:,7], label="ATP", color="gray")
# ax6.legend()
# ax6.set_ylim(0, 1.5)
# ax7 = plt.subplot(717)
# ax7.plot(data[:,0], data[:,7], label="ATP", color="gray")
# ax7.plot(data[:,0], data[:,10], label="PKaiC up", color="blue")
# ax7.plot(data[:,0], data[:,11], label="PKaiC dn", color="green")
# ax7.legend()
# ax7.set_xlabel("Time (h)")

fig.tight_layout()
filename = "Limit Cycle, eta = " + str(eta) + ", ATP = " + str((ATP1, ATP2)) + ", stochastic = " + str(stoch_flag) + ", env_noise = " + str(env_noise)
fig.savefig("Signal and P-KaiC traces, " + filename + ".pdf")
# fig.savefig('test.pdf')
print("Program run complete.")

print("Runtime:", endTime - startTime, "seconds.")


plt.figure()
plt.plot(data[initial_index_of_data:,10], data[initial_index_of_data:,11])
plt.xlabel('PKaiC up')
plt.ylabel('PKaiC dn')
plt.savefig(filename + ".pdf")

time = data[:,0]
x    = data[:,10]
y    = data[:,11]


t_interpolate = np.arange(starttime, endtime, dt)
y_interpolate = ladder_interpolate(time, y, t_interpolate)
x_interpolate = ladder_interpolate(time, x, t_interpolate)

time = data[initial_index_of_data:,0]
x    = data[initial_index_of_data:,10]
y    = data[initial_index_of_data:,11]

# old_answer = find_mi(time, x, y)
# answer = find_mi(t_interpolate, x_interpolate, y_interpolate)
# print('MI = ' + str(answer) + ', old MI = ' + str(old_answer))

answer = find_mi(t_interpolate, x_interpolate, y_interpolate)
print('MI = ' + str(answer))





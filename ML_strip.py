from Data_Importing import import_data
from EC_MS import numerize, synchronize, time_cut, plot_masses
# from Isotope_Ratio import predict_M34
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy
import Chem  #physical constants and Molar Masses
from scipy.integrate import odeint


def get_monolayer_signal_area(CA_and_MS, N_sites = 4e13, verbose = 1, tspan = 0 ):
    '''
    This will calculate the area expected under a QMS signal from desorption of 
    one monolayer of oxygen based on the M32 signal and associated electrical 
    current at steady-state oxidation of normal water and information about the 
    catalytic surface (simplest is the number of sites)
    '''
    if verbose:
        print('\n\nfunction \'get_monolayer_signal_area\' at your command!')
        
    #cut off the first half of the data to deal with steady state only
    if not tspan:
        tspan = CA_and_MS['tspan_2']
        tspan[0] = (tspan[1] + tspan[0])/2
    Data = CA_and_MS.copy()
        #deepcopy needed to not fuck up CA_and_MS elsewhere
    Data = time_cut(Data, tspan)    

    
    if verbose:
        print('plotting steady state')
        plot_masses(Data, Colors = {'M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m'})
    #get the useful EC current and MS signal variables
    I_x = Data['time/s']           
    I_y = Data['I/mA']             # in mA
    s_x = Data['M32-x']
    s_y = Data['M32-y'] + Data['M34-y'] + Data['M36-y']
    
    #integrate the current and the signal
    Q_CA = np.trapz(I_x,I_y)*1e-3       #charge passed during second half of CA, in C
    A_CA = np.trapz(s_x,s_y)            #integrated QMS signal (also in C, incidentally)
    
    #charge expected from one turnover per site:
    Q_ML = N_sites*4*Chem.qe            #charge of 4 electrons per surface site, in C
    
    #signal area expected from one monolayer
    A_ML = Q_ML*A_CA/Q_CA               #expected area of one monolayer
    
    if verbose:
        print(str(Q_CA) + ' C passed for an integrated M32 signal of ' + str(A_CA) + '\n' + str(Q_ML) + ' C in a monolayer, corresponds to M32 signal of ' + str(A_ML))
    return A_ML


def diffusion_ode(c_vec, t, alpha, pulse_fun):
    '''
    sniffer ode problem. returns rate of change of c at points throught working 
    volume (dc_dt) ad function of their present values (c_vec), the system 
    parameter alpha = L*h/D, and the flux from the electrode
    (pulse_fun(t))
    
    Non-dimensionalized, as described Scott's thesis, chapter 2
    '''
    N = len(c_vec)      #level of discretization
    dx = 1/N            #non-dimensional interpoint distance
    

    c_0 = c_vec[0]*(1 - alpha*dx)                     
    #boundary condition at membrane, dc/dx = alpha*c0
    c_end = c_vec[-1] + pulse_fun(t)*dx          
    #boundary condition at electrode, dc/dx = - j_end
    #size of pulse_fun is arbitrary. will normalize after solving.    
        
    c_vec_minus = np.append(c_0, c_vec[:-1])
    c_vec_plus = np.append(c_vec[1:], c_end)
    
    dc_dt = (c_vec_plus - 2*c_vec + c_vec_minus)/dx**2 #rate of change due to diffusion
    return dc_dt
    
 
def solve_diffusion_ode(ax1 = 'new', Area = 9.93e-10, t_pulse = 0.001, verbose = 1, 
    spec = 'k--', plot_label = 'predicted ML signal', logplot = 1,
    N = 30 ,             #discretization in x
    L = 100e-6  ,        #distance between electrode and membrane in m
    D_O2 = 2.10e-9 ,     #diffusion constant of oxygen in water in m^2/s
    h = 6.4e-5          #mass transport coefficient at sniffer membrane, from Master's thesis, m/s 
    ):
    '''
    Solves sniffer diffusion ode for a pulse as described in Scott's master 
    thesis chapter 2.
    (First time I solve a partial differential equation in python)
    '''
    if verbose:
        print('\n\nfunction \'solve_diffusion_ode\' at your command! \nGenerating ' + plot_label)
    alpha = L*h/D_O2
    if verbose:
        print('\tsystem parameter alpha = ' + str(alpha))
    c_0 = np.zeros([N])
    t0 = L**2/D_O2    
    if verbose:
        print('\tnon-dimensionalizing on time with \n\t t_0 = ' + str(t0) + ' s')
    tau_pulse = t_pulse/t0   #non-dimensionalized pulse length
    tau_vec = np.linspace(0, 10 + tau_pulse, 300)
    
    def pulse_fun(tau):       #two-sided step function simulating an electrochemical pulse
        if tau>0 and tau<tau_pulse:
            return 1
        return 0
    
    ##solve the ode!
    Cc = odeint(diffusion_ode, c_0, tau_vec, args=(alpha, pulse_fun) )    
    
    s = Cc[:,0]                 #signal is proportional to concentration at membrane
    t = t0*tau_vec  
    A_out = np.trapz(s,t)
    A_in = tau_pulse*N/L  #to test the diffusion problem
    if verbose:
        print('A_in = ' + str(A_in) + '\nA_out = ' + str(A_out) + '\nnormalizing to A_ML = ' + str(Area))
    s = Area/A_out*s
    


    if verbose:
        if ax1 == 'new':
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.set_xlabel('t / s')
            y_string = 'singal / [A]'
            if logplot:
                ax1.set_ylabel('log(' + y_string + ')')
        else:
            ylim = ax1.get_ylim()
            if logplot:
                s = s + np.exp(np.log(10)*-12.3) 
            #rudimentary background addition, so that I don't get problems on a log scale
        if logplot:
            s = np.log(s)/np.log(10)
        ax1.plot(t,s, spec, label = plot_label)
        print('just plotted, I think.')
        ax1.legend()
    
    return [t,s]

if __name__ == '__main__':
    
    #plt.close()    
    
    default_directory = '/home/soren/Desktop/Sniffer_Experiments/O18_NiNPs/00_python/test_files/'    

    MS_file = default_directory + 'QMS_data.txt'
    CA_file = default_directory + '01_O16_10_CA_C01.mpt'
    #CA_file = default_directory + '02_O16_to_O18_10_CA_C01.mpt'
    #CA_file = default_directory + '03_O18_10_CA_C01.mpt'
    #CA_file = default_directory + '04_O18_to_O16_10_CA_C01.mpt'


    MS_Data = numerize(import_data(MS_file,data_type = 'MS'))
    CA_Data = numerize(import_data(CA_file))
    
    CA_and_MS = synchronize([MS_Data, CA_Data])
    
    tspan0 = CA_and_MS['tspan_2']    #'tspan' is given in seconds since midnight
                                    #'tspan_2' is given in seconds since start of CA file
                                    #'tspan_1' is given in seconds since start of QMS file                             
    tspan = [0,0]
    tspan[1] = tspan0[1] + 100 + 0*(tspan0[1]-tspan0[0])
    tspan[0] = tspan0[0] - 20 - 0*(tspan0[1]-tspan0[0]) 
    
    Area = get_monolayer_signal_area(CA_and_MS, verbose = 0, 
                                     #tspan = [1000, 1500]) # for O18_to_O16
                                     )
    CA_and_MS = time_cut(CA_and_MS,tspan)
    
    plot_masses(CA_and_MS, logplot = 1,
     #           Colors = {'M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m', 'predicted M34':'--k'})
                 Colors = {'M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m'})
    #MS_Data_1 = predict_M34(MS_Data_1)  
    

    solve_diffusion_ode(ax1 = plt.gca(), Area=Area, t_pulse = 1, logplot = 1)
    
  
    
    
    
    
    
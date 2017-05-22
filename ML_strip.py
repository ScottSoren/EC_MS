

"""
Most recently edited: 16I23

@author: Scott

This module has functions specifically designed to investigate 
single-turnover phenomena in EC-MS.
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import os

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    from Plotting import plot_masses
    import Chem
else:                           #then we use relative import
    from .Plotting import plot_masses
    from . import Chem


def get_monolayer_signal_area(CA_and_MS, N_sites = 4.0e13, verbose = 1, 
                              tspan = 0, masses = ['M32','M34','M36'],
                              Colors = {'M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m'},  
                              n_el = 4, n_el_test = 4,
                              FA = 1, sensitivity_ratio = 1):
    '''
    This will calculate the area expected under a QMS signal from 
    'n_el' -electron desorption of one monolayer of a product based on:
      (1) the signal and associated electrical current in 'CA_and_MS' during 'tspan', 
          under steady-state 'n_el_test' -electron production with
          faradaic efficiency 'FA' of a test product specified in 'masses';
      (2) 'sensitivity_ratio' between the product and the test product;
      (3) information about the catalytic surface where the
          simplest is the total number of sites 'N_sites'
    Defaults are for OER isotope experiments
    '''
    if verbose:
        print('\n\nfunction \'get_monolayer_signal_area\' at your service!')

        
    #cut off the first half of the data to deal with steady state only
    if not tspan:
        tspan = CA_and_MS['tspan_2'].copy()
        tspan[0] = (tspan[1] + tspan[0])/2
    Data = CA_and_MS.copy()
        #deepcopy needed to not fuck up CA_and_MS elsewhere
    Data = time_cut(Data, tspan)    
   
    if verbose:
        print('plotting steady state')
        ax1 = plot_masses(Data, Colors = Colors)
        ax1.set_title('Steady state used to get monolayer signal area')
    #get the useful EC current and MS signal variables
    I_x = Data['time/s']           
    I_y = Data['I/mA']             # in mA
    if type(masses)==list:
        s_x = Data[masses[0]+'-x']
        s_y = np.zeros(np.shape(s_x))
        for mass in masses:
            s_y = s_y + Data[mass + '-y']
    else:
        s_x = Data[masses +'-x']
        s_y = Data[masses +'-y']
    
    #integrate the current and the signal
    Q_CA = np.trapz(I_y,I_x)*1e-3*FA       #charge passed to test product during second half of CA, in C
    A_CA = np.trapz(s_y,s_x)               #integrated QMS signal (also in C, incidentally)
    
    # number of molecules of test product produced during CA
    N_CA = Q_CA/(n_el_test * Chem.qe)
    
    #signal area expected from one monolayer
    A_ML = A_CA *N_sites/N_CA* sensitivity_ratio    #expected area of one monolayer
    
    if verbose:
        print('Integrated signal of ' + str(A_CA) + ' observed for passage of ' + str(Q_CA) +
            ' Coulumbs, corresponding to ' + str(N_CA) + ' molecules of test product.\n' +
            'Integrated product signal area expected for desorption of 1 ML = ' + str(N_sites) + ' sites \n' +
            '  is: ' + str(A_ML))
        print('function \'get_monolayer_signal_area\' finished!\n\n')    
    
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
    #magnitude of pulse_fun is arbitrary. will normalize after solving.    
        
    c_vec_minus = np.append(c_0, c_vec[:-1])
    c_vec_plus = np.append(c_vec[1:], c_end)
    
    dc_dt = (c_vec_plus - 2*c_vec + c_vec_minus)/dx**2 #rate of change due to diffusion
    
    return dc_dt
    
 
def solve_diffusion_ode(ax1 = 'new', Area = 9.93e-10, t_pulse = 0.001, verbose = 1, 
    spec = 'k--', plot_label = 'predicted ML signal', logplot = 1, background = 6.0e-13,
    N = 30 ,             #discretization in x
    L = 100e-6  ,        #distance between electrode and membrane in m
    D = 2.10e-9 ,        #diffusion constant of oxygen in water in m^2/s
    h = 6.4e-5           #mass transport coefficient at sniffer membrane, from Master's thesis, m/s 
    ):
    '''
    Solves sniffer diffusion ode for a pulse as described in Scott's master 
    thesis chapter 2.
    (First time I solve a partial differential equation in python)
    '''
    if verbose:
        print('\n\nfunction \'solve_diffusion_ode\' at your service! \nGenerating ' + plot_label)
    alpha = L*h/D
    if verbose:
        print('\tsystem parameter alpha = ' + str(alpha))
    c_0 = np.zeros([N])
    t0 = L**2/D    
    if verbose:
        print('\tnon-dimensionalizing on time with \n\t t_0 = ' + str(t0) + ' s')
    tau_pulse = t_pulse/t0   #non-dimensionalized pulse length
    tau_vec = np.linspace(0, 5 + tau_pulse, 300)
    
    def pulse_fun(tau):       #two-sided step function simulating an electrochemical pulse
        if tau>0 and tau<tau_pulse:
            return 1
        return 0
    
    ##solve the ode!
    Cc = odeint(diffusion_ode, c_0, tau_vec, args=(alpha, pulse_fun) )    
    
    s = Cc[:,0]                 #signal is proportional to concentration at membrane
    t = t0*tau_vec  
    A_out = np.trapz(s,t)
    A_in = tau_pulse            #to check the diffusion problem, 
    #but I can't actually figure out what this should be given the poor way I've non-dimensionalized this.
    if verbose:
        print('A_in = ' + str(A_in) + '\nA_out = ' + str(A_out) + '\nnormalizing to A_ML = ' + str(Area))
    s = Area/A_out*s
    
    s = s + background     # add the QMS intrinsic background to signal

    if verbose:
        if ax1 == 'new':
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.set_xlabel('t / s')
            y_string = 'singal / [A]'
            ax1.set_ylabel(y_string)
            xlim = 0
        else:
            xlim = ax1.get_xlim()

        ax1.plot(t,s, spec, label = plot_label)
        print('just plotted, I think.')
#        ax1.legend()
        if logplot:
            ax1.set_yscale('log')
        if xlim:
            ax1.set_xlim(xlim)
    if verbose:
        print('function \'solve_diffusion_ode\' finished!')
    return [t,s]




if __name__ == '__main__':
    
    from Data_Importing import import_data
    from Combining import numerize, synchronize, time_cut, plot_masses, plot_masses_and_I
    from EC import select_cycles
    #plt.close()    
    
    default_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))  

    #MS_Data = numerize(import_data(MS_file,data_type = 'MS'))
    import_raw_data = 0
    if import_raw_data:
        MS_Data_0 = import_data(default_directory, data_type='MS')
        CA_Data_0 = import_data(default_directory, data_type='EC')

    cycles = np.array([17])
    
    #get the time span of the pulse from select cycles
    tspan0 = select_cycles(CA_Data, cycles)['tspan'] 
    
    #synchronize the entire files, with t=0 at the start of the pulse
    CA_and_MS = synchronize([MS_Data, CA_Data], cutit = 0, t_zero = tspan0[0])    
    
    #define a span from a bit before to a bit after the pulse, which starts at t=0
    tspan = [0,0]
    tspan[0] =  -20
    tspan[1] = tspan0[1] - tspan0[0] + 50
    #cut the data accordingly
    CA_and_MS = time_cut(deepcopy(CA_and_MS),tspan)      
    
    n_el_test = -2    #number of electrons in HER
    n_el = -6          #number of electrons for CORR to methane
    sensitivity_ratio = 3.419/0.969    #ratio of sensitivity factors on M15 for methane to M2 for H2
    # from NIST: http://physics.nist.gov/PhysRefData/Ionization/molTable.html
    
    a0_Cu = 3.61e-10 #lattice constant for Cu / m
    site_density = 4/(np.sqrt(3)*a0_Cu**2) #density of surface sites on Cu(111) / m^-2
    A_electrode = 0.2e-4    #electrode area, m^2
    coverage = 2*0.05       #SA_Cu/SA_electrode at 5% projected coverage and assuming hemispherical NPs
    N_sites = A_electrode*coverage*site_density    

    tspan_SS = [0.5*(tspan0[1]-tspan0[0]), tspan0[1]-tspan0[0]]    
    #steady state tspan to calibrate with M2 signal during CA / s
    
    
    Colors = {'M2':'b', 'M15':'r', 'M27':'g'}
    A_ML = get_monolayer_signal_area(CA_and_MS, N_sites = N_sites, verbose = 0, 
                              tspan = tspan_SS, masses = ['M2'], Colors = Colors,
                              n_el = n_el, n_el_test = n_el_test,
                              FA = 1, sensitivity_ratio = sensitivity_ratio)
    
    D_CH4 =  1.49e-9        #diffusion constant of methane in water / [m^2/s]
    
    #and plot it!
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(121)
    plot_masses(CA_and_MS, 
                Colors = Colors,
     #          Colors = {'M2':'k','M32':'b', 'M34':'r', 'M36':'g', 'M18':'c', 'M20':'m'})
                ax1 = ax1
                )
                
    [t,s] = solve_diffusion_ode(Area = A_ML, t_pulse = 5, verbose = 0, 
                        spec = 'k--', plot_label = 'predicted ML signal', logplot = 1, 
                        background = 2.0e-12,
                        N = 30 ,             #discretization in x
                        L = 100e-6  ,        #distance between electrode and membrane in m
                        D = D_CH4 ,        #diffusion constant of oxygen in water in m^2/s
                        h = 3e-5           #mass transport coefficient at sniffer membrane, from Master's thesis, m/s 
                        )
    t = t + 3 #apply offset to account for the shitty slow QMS
    
    ax1.plot(t,s,'k--',label = 'predicted')        
    ax1.legend_.remove()
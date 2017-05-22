# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:50:10 2016
Most recently edited: 16J27

@author: scott

This module will include functions for modelling mass transport
and fitting time response in the EC-MS setup. 

See Scott's MSc thesis, chapter 2 and section 3.3, for discussion inc. prior
implementation in Matlab


"""


from __future__ import division, print_function

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import odeint

if os.path.split(os.getcwd())[1] == 'EC_MS':      
                                #then we're running from inside the package
    import Chem
    from Molecules import Molecule
else:                           #then we use relative import
    from . import Chem
    from .Molecules import Molecule


def fit_exponential(t,y):
    '''
    A previous attempt at this had used scipy.optimize.minimize.
    '''
    t = t - t[0]        #zero time axis
    tau_i = t[-1]/10      #guess at time constant 
    #tau_i = t[-1]      #often can't solve with this guess. A smaller tau helps.
    y0_i = y[-1]        #guess at approach value
    y1_i = y[0]         #guess at true initial value
    pars_i = [tau_i, y0_i, y1_i]    

    def exp_fun(x, tau, y0, y1):
        z = y0 + (y1 - y0) * np.exp(-x / tau)
#        print([tau,y0,y1]) #for dianosin curve_fit problems
        return z

    pars, pcov = curve_fit(exp_fun, t, y, p0=pars_i)
#    pars = [tau, y0, y1]

    return pars

    
def fit_step(t, y, tpulse=0, step='fall', 
             ax=None, spec='r--', label=None, verbose=1):
    '''
    use step='rise' to fit onset and step='fall' to fit tail
    assumes that data starts from the start of the pulse
    '''
    if verbose:
        print('\n\nfunction \'fit_step\' at your service!\n')
    
    #zero time axis
    t0 = t[0]
    t = t - t0
    print('t0 = ' + str(t0))
    if type(t) is list:
        t = np.array(t)         #17B02
    if step=='fall':   
        I_tail = np.array([I for (I,t_I) in enumerate(t) if tpulse<t_I])
 #       print(I_tail)
        t_fit = t[I_tail] - tpulse
    elif step =='rise':
        if tpulse == 0:
            tpulse = t[-1]
        I_tail = np.array([I for (I,t_I) in enumerate(t) if t_I<tpulse])
        t_fit = t[I_tail]
    else:
        print('use step=\'rise\' to fit onset and step=\'fall\' to fit tail')

    pars = fit_exponential(t_fit, y[I_tail])    
    
    if ax:
        tau = pars[0]
        y0 = pars[1]
        y1 = pars[2]    
        y_fit = y0 + (y1 - y0)*np.exp(-t_fit/tau)
        t_fit = t_fit + t0          #put time axis back
        if step =='fall':
            t_fit = t_fit + tpulse
        plot1 = ax.plot(t_fit, y_fit, spec, label=label)
        if label:
            if label == 'tau':
                label = 'tau = {0:5.2f} s'.format(tau)
            I_text = int(len(t_fit)/2)
            t_text = t_fit[I_text]
            y_text = y_fit[I_text]
            ax.text(t_text, y_text, label, color=spec[0])
        
        
    if verbose:
        print('tau = ' + str(pars[0]) + ' s')
        print('\nfunction \'fit_step\' finished!\n\n')
    return pars
    

    
def stagnant_diffusion_ode(C, T, pars):  #Note that C comes before T here!
    '''
    returns rate of change dC/dT of concentration profile for 
    non-dimensionalized stagnant sniffer diffusion problem.
    C = C(X) where X goes from 0 (membrane) to 1 (electrode)
    T is time non-dimensionalized on the diffusion scale t0 = L²/D
    pars:
    [0] alpha = h*L/D is the system parameter.
    [1] J_fun returns the flux from the electrode as a unction of T. The flux 
    scale J0 is used to define the concentration scale, C0 = J0*L/D
    
    #modified 17B02 to enable carrier gas introduction of element using Cg
    '''
    alpha = pars[0]
    J_fun = pars[1]
    Cg = pars[2]
    N = np.size(C)
    dX = 1 / N    #N+1? N-1? I can never remember what's most 'correct'
    C_N = C[-1] + J_fun(T) * dX      # boundary condition dC/dX = J(T) at electrode
    C_ = C[0] - alpha * (C[0] - Cg) *dX       # boundary condition dC/dX = -alpha*(C-Cg) at membrane
    C_up = np.append(C[1:], C_N)
    C_down = np.append(C_, C[:-1])
    d2CdX2 = (C_up - 2*C + C_down)/(dX*dX)  #second derivative of C wrt X
    dCdT = d2CdX2       #Fick's second law
    return dCdT    

def solve_stagnant(pars, Tspan, startstate, N=30, flux=1, verbose=1):
    '''solves the stagnant sniffer partial differential equations. 
    pars[0][0] is alpha = h*L/D is the system parameter.
    pars[0][1] is J_fun. Returns the flux from the electrode as a unction of T. 
    Tspan is [Tstart, Tfinish] on the diffusion timescale t0 = L²/D
    C0 = C0(X) is start concentration profile. If size(C0) = 1, assumes a 
    uniform concentration.
    N is descretization (read from C0)
    flux = 0 to return entire concentration profile (on C0 = J0*L/D scale)
    flux = 1 to return flux through membrane (on J0 scale)
    '''
    if verbose:
        print('\n\nfunction \'solve_stagnant\' at your service!\n')
    
    if startstate == 'zero':
        C0 = np.zeros([N])
    elif startstate == 'steady':
        alpha = pars[0][0] #for some reason pars is doubly packed.
        C0 = 1/alpha + np.linspace(0,N)/(N+1) #Scott's MSc thesis, p. 53
    elif startstate == 'saturated':
        Cg = pars[0][2]
        C0 = np.ones([N]) * Cg
    elif np.size(startstate) == 1:
        C0 = np.ones([N]) * startstate
    else:
        C0 = startstate
        N = np.size()   
        
    if np.size(Tspan) == 2:
        Tspan = np.linspace(Tspan[0], Tspan[1], 100)
    CC = odeint(stagnant_diffusion_ode, C0, Tspan, args=pars)    
    #16J18_02h10: this crashes the kernel, I don't know why... 18h56 found it! c before t!
    
    J = (CC[:,1] - CC[:,0]) * N  # (positive) J = dC/dX with dC = C0 - C_ and dX = 1 / N 
     
    if verbose:
        print('solution shape: ' + str(np.shape(CC)))
        print('\nfunction \'solve_stagnant\' finished!\n\n')    
    if flux:
        return Tspan, J
    else:
        return Tspan, CC
          
        
def stagnant_pulse(tj=None, tpulse=10, tspan=None, j_el = -5,
                   L=100e-6, A=0.196e-4, q0=1e15/Chem.NA, p_m=1e5,
                   mol='H2', p_gas=0, normalize=False, 
                   D=None, kH=None, n_el=None, Temp=None, 
                   unit = 'umol / (m^2*s)', flux_direction='out',
                   verbose=1, plot_type='flux', startstate='zero'):
    '''                   
    Models a pulse of current towards a specified product in our EC-MS setup.
    Theory in chapter 2 of Scott's masters thesis.
    all arguments are in pure SI units. The electrode output can either be given
    as a steady-state square pulse of electrical current (tpulse, j_el, n_el), 
    or as a measured current (tj[1]) as a function of time (tj[0])
    
    #17B02: p_gas is the partial pressure of the analyte in the carrier gas.
    # this enables, e.g., CO depletion modelling.
    '''    
    if verbose:
        print('\n\nfunction \'stagnant_pulse\' at your service!\n') 
    if type(mol) is str:            
        mol = Molecule(mol)
    if Temp is not None:
        mol.set_temperature(Temp)
    else:
        Temp = 298.15       #standard temperature in K
    if D is None:
        D = mol.D
    if kH is None:
        kH = mol.kH
    if n_el is None and not normalize:
        n_el = mol.n_el                 
    if tspan is None:
        if tj is None:
            tspan = [-0.1*tpulse, 1.2*tpulse]
        else:
            tspan = [tj[0][0], tj[0][-1]]

    h = kH*Chem.R*Temp*q0/(p_m*A)  #mass transfer coefficeint

    alpha = L*h/D               #system parameter

    #non-dimensional scales:
    t0 = L**2/D
    if tj is None:
        if normalize:
            j0 = 1
        else:
            j0 = j_el / (n_el * Chem.Far)
    else:
        t = tj[0]
        if normalize:
            j0 = 1
            j = tj[1]/np.max(np.abs(tj[1]))
        else:
            j = tj[1] / (n_el * Chem.Far) # A/m^2 --> mol/(m^2*s) 
            j0 = max(np.abs(j)) 
    c0 = j0*L/D
    tau = L**2/(2*D) + L/h       
    #from the approximate analytical solution, Scott's thesis appendix D
    
    Tpulse = tpulse/t0
    Tspan = np.linspace(tspan[0],tspan[1],1000)/t0
    
    if tj is None: 
        def J_fun(T):
            if T < 0: 
                return 0
            if T < Tpulse:
                return 1
            return 0
    else:
        T_in = t / t0
        J_in = j / max(np.abs(j))
        print('max(J_in) = ' + str(max(J_in)))
        def J_fun(T):
            if T < T_in[0]:       #assume no current outside of the input tj data
                return 0
            if T < T_in[-1]:
                return np.interp(T, T_in, J_in) 
            return 0
    
    c_gas = p_gas / (Chem.R*Temp)
    cg = c_gas / kH  #concentration analyte in equilibrium with carrier gas, 17B02
    Cg = cg / c0 #non-dimensionalized concentration analyte at equilibrium with carrier gas, 17B02
    
    pars = ([alpha, J_fun, Cg],) #odeint needs the tuple. 17A12: Why ?!
    
    [T, CC] = solve_stagnant(pars, Tspan, startstate, flux=False)   
    cc = CC*c0
    t = T*t0
    j = h * (cc[:,0] - cg)       #mass transport at the membrane
    #j1 = D * (cc[:,1] - cc[:,0])        
        #fick's first law at the membrane gives the same j :)
    if verbose:
        print('q0 = ' + str(q0) + ' mol/s, h = ' + str(h) + ' m/s, alpha = ' + str(alpha) + 
            ', j0 = ' + str(j0) + ' mol/(m^2*s), max(j)/j0 = ' + str(max(j)/j0) + 
            ',  t0 = ' + str(t0) + ' s, c0 = ' + str(c0) + ' mol/m^3' + 
            ', tau (analytical) = ' + str(tau) + ' s' + ', cg = ' + str(cg) + ' mM')  
        
    # get ready to plot:
    N = np.shape(cc)[1]
    x = np.arange(N)/(N-1)*L          
        #this will only be used for heatmap, so it's okay if dx isn't quite right.
    if 'cm^2' not in unit:
        j = j*A
    if unit[0] == 'u':
        j = j*1e6
    elif unit[0] == 'n':
        j = j*1e9
    elif unit[0] == 'p':
        j = j*1e12
    if flux_direction=='in':
        j = -j
    axes = None
    
    # and plot! 
    if plot_type == 'flux':    
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.plot(t, j, label='simulated flux')
        ax1.set_xlabel('time / s')
        ax1.set_ylabel('flux / [' + unit + ']')
        axes = ax1
        
    elif plot_type == 'heat' or plot_type == 'both':
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        
        #t_mesh, x_mesh = np.meshgrid(t,x)
        #img = ax1.contourf(t_mesh, x_mesh*1e6, np.transpose(cc,[1,0]), cmap='Spectral', 
        #                   levels=np.linspace(np.min(cc),np.max(cc),100))
        
        # imshow objects seem more versatile than contourf for some reason.
        
        xrange = [min(t), max(t)]
        yrange = [min(x*1e6), max(x*1e6)]
        img = ax1.imshow(np.transpose(cc,[1,0]), 
                         extent=xrange[:] + yrange[:],  #have to be lists here!
                        aspect='auto', origin='lower',
                         cmap = 'Spectral')        
        
        cbar = plt.colorbar(img, ax=ax1)
        cbar.set_label('concentration / mM')
        ax1.set_xlabel('time / s')
        ax1.set_ylabel('position  / um')
        
#        print('plot_type = ' + plot_type)
        if plot_type == 'both':
            ax2 = ax1.twinx()
            ax2.set_ylabel('flux / [' + unit + ']')
            ax2.plot(t, j, 'k-')
            cbar.remove()
            ax3 = img.figure.add_axes([0.85, 0.1, 0.03, 0.8])
            cbar = plt.colorbar(img, cax=ax3)
            cbar.set_label('concentration / mM')
            ax1.set_xlim(tspan)
            print('returning three axes!')
            axes = [ax1, ax2, ax3]
        else:
            axes = [ax1, cbar]
        
    if normalize:
        s_int = np.trapz(j,t)
        if verbose:
            print('normalizing from area = ' + str(s_int))
        j = j / s_int
    
    if verbose:  
        print('\nfunction \'stagnant_pulse\' finished!\n\n')  
    
    if axes is not None:
        return (t, j, axes)
    return (t, j)
   

def delta_response(L=100e-6, q0=1e15/Chem.NA, 
                   mol='H2', D=None, kH=None, n_el=None, 
                   A=0.196e-4, p_m=1e5, Temp=298.15,
                   verbose=True, tspan='auto', N_t=1000):
    '''
    Returns the normalized response of a delta-function input as (t, j).
    There's probably a much smarter way to do it, but for now I'll just do
    a millisecond pulse
    '''
    if D is None or kH is None:
        if type(mol) is str:
            mol = Molecule(mol)
        if D is None:
            D = mol.D
        if kH is None:
            kH = mol.kH
    if verbose:
        print('calculating a delta function response.')
    h = kH*Chem.R*Temp*q0/(p_m*A)  #mass transfer coefficeint
    tau = L/h + L**2/(2*D)
    if tspan=='auto':
        tspan = [0, 4*tau]
    
    t = np.linspace(tspan[0], tspan[1], N_t)
    j = np.append(np.array([1]), np.zeros(N_t-1))
    tj = [t,j]
    print(type(tj))
    return stagnant_pulse(tj=tj, normalize=True, tspan=tspan,
                   L=L, A=A, q0=q0, p_m=p_m,
                   D=D, kH=kH, n_el=n_el, Temp=Temp, 
                   verbose=True, plot_type=None)
    
    
                   
   
if __name__ == '__main__':

    
    from Data_Importing import import_data
    from Combining import synchronize
    from EC import select_cycles, plot_CV_cycles, plot_vs_time
    from Plotting import plot_masses_and_I
  

    plt.close('all')  
    
    ## Test and fit the numerical modelling of pulse

    tpulse = 10
    
    t, j, *ax = stagnant_pulse(tpulse=tpulse, tspan=[0, tpulse*2], j_el=5, 
                          plot_type='both',
                          mol='O2')        
    ax2 = ax[1]
    pars1 = fit_step(t, j*1e6, tpulse=tpulse, ax=ax2, step='rise')
    pars2 = fit_step(t, j*1e6, tpulse=tpulse, ax=ax2, step='fall')



        
        
        ##
    
    
    
    




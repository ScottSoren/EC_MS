# -*- coding: utf-8 -*-
"""
Created on Wed April 22 19:50:34 2020
Most recently edited: 20D22

@author: Kevin

This module contains the mass_transport class that allows to back calculate
(deconvolute) reaction rates at the electrode-electrolyte interface and vice versa.
EC_to_MS is equivalent to the stagnant_operator function in Time_Response.py but
uses convolution and a numerical implementation of the inverse Laplace transform
to solve the problem instead of an ODE.
MS_to_EC calculates reaction rates based on a MS signal and is the core function
of this module. The most reliable results are obtained when "calibrating" the
mass transport by fitting of impulse responses of the MS-signal implemented in
fit_impulse.
The term kernel is used interchangebly with the term impulse response.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import signal
from mpmath import invertlaplace, sinh, cosh, sqrt, exp, erfc, pi, tanh, coth
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft, ifftshift, fftfreq


from .Molecules import Molecule
from . import Chem




def gen_kernel(x, D, L, h, c, fourier = False, matrix = False, norm=False, verbose = False):
    """
    Returns the kernel along x (time or frequency) with
        D: Diffusion coefficient
        L: Distance between WE and chip
        h: mass transport coefficient at liquid/gas interface
        c: pumping speed correction parameter
    In time domain the output can also be in form of the convolution matrix generated
    from the kernel.
    """
    if verbose:
        print('D: ' + str(D) + ' || L: ' + str(L) + ' || h: ' + str(h) + ' || c: ' + str(c))

    if fourier:
        print('Kernel in fourier space implemented soon')
        fs = lambda s: 1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))#*1/(2*c*sqrt(2*pi()))*exp(sqrt(s/c)/8)*erfc(s/(2*sqrt(2)*c))

    else:
        tdiff=x*D/(L**2)
        a=L*h/D
        fs = lambda s: 1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))*1/(s+c*L**2/D)
        #0.5*sqrt(pi()/c)*(1-erf(s/(2*sqrt(pi()))))*exp(s**2/(4*c))
        #*1/(s+c)
        #1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))*

        kernel=np.zeros(len(x))
        for i in range(len(x)):
            kernel[i]=invertlaplace(fs,tdiff[i],method='talbot')


        #disp = taylor_dispersion(x, D=c, L=l_cap)

        #kernel = signal.convolve(kernel, disp, mode='full')
        #kernel = kernel[:-(len(x)-1)]

        if matrix:
            kernel = np.tile(kernel,(len(kernel),1))
            i=1
            while i < len(x):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i=i+1

    if norm:
        area = np.trapz(kernel, x)
        return abs(kernel/area)
    else:
        return kernel

def gen_kernel2(x, D, L, V, V_dot, kH, fourier = False, matrix = False, norm=False, verbose = False):
    """
    Returns the kernel along x (time or frequency) with
        D: Diffusion coefficient
        L: Distance between WE and chip
        h: mass transport coefficient at liquid/gas interface
        c: pumping speed correction parameter
    In time domain the output can also be in form of the convolution matrix generated
    from the kernel.
    """
    if verbose:
        print('D: ' + str(D) + ' || L: ' + str(L) + ' || h: ' + str(h) + ' || c: ' + str(c))

    if fourier:
        print('Kernel in fourier space implemented soon')
        fs = lambda s: 1/(sqrt(s)*sinh(sqrt(s))+(V*kH/(0.196e-4)/L)*(s+V_dot/V*L**2/D)*cosh(sqrt(s)))#*1/(2*c*sqrt(2*pi()))*exp(sqrt(s/c)/8)*erfc(s/(2*sqrt(2)*c))

    else:
        tdiff=x*D/(L**2)




        fs = lambda s: 1/(sqrt(s)*sinh(sqrt(s))+(V*kH/(0.196e-4)/L)*(s+V_dot/V*L**2/D)*cosh(sqrt(s)))
        #0.5*sqrt(pi()/c)*(1-erf(s/(2*sqrt(pi()))))*exp(s**2/(4*c))
        #*1/(s+c)
        #1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))*

        kernel=np.zeros(len(x))
        for i in range(len(x)):
            kernel[i]=invertlaplace(fs,tdiff[i],method='talbot')


        #disp = taylor_dispersion(x, D=c, L=l_cap)

        #kernel = signal.convolve(kernel, disp, mode='full')
        #kernel = kernel[:-(len(x)-1)]

        if matrix:
            kernel = np.tile(kernel,(len(kernel),1))
            i=1
            while i < len(x):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i=i+1

    if norm:
        area = np.trapz(kernel, x)
        return abs(kernel/area)
    else:
        return kernel

def gen_kernel3(x, D_l, L, V_dot, D_g, l, kH, fourier = False, matrix = False, norm=False, verbose = False):
    """
    Returns the kernel along x (time or frequency) with
        D: Diffusion coefficient
        L: Distance between WE and chip
        h: mass transport coefficient at liquid/gas interface
        c: pumping speed correction parameter
    In time domain the output can also be in form of the convolution matrix generated
    from the kernel.
    """
    if verbose:
        print('D: ' + str(D) + ' || L: ' + str(L) + ' || h: ' + str(h) + ' || c: ' + str(c))

    if fourier:
        print('Kernel in fourier space implemented soon')
        fs = lambda s: 1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))#*1/(2*c*sqrt(2*pi()))*exp(sqrt(s/c)/8)*erfc(s/(2*sqrt(2)*c))

    else:
        #tdiff=x*D_l/(L**2)
        #alpha = V_dot*l_diff/(0.196e-4*D_g)
        #beta= D_g/D_l*L/l_diff
        #gamma= D_l/D_g*l_diff**2/L**2
        V = V_dot/(0.196e-4)
        eps=0.05

        fs = (
            lambda s: -4*kH*sqrt(D_g*s)*sqrt(D_l*s)*exp(L*sqrt(s/D_l) + l*sqrt(s/D_g))/(

                (V - sqrt(D_g*s))
                *(
                    2*D_l*s*exp(2*L*sqrt(s/D_l))
                    - sqrt(D_l*s)*(kH*sqrt(D_g*s)
                    + sqrt(D_l*s))*(exp(2*L*sqrt(s/D_l)) + 1)
                )
                -(V + sqrt(D_g*s))*(2*D_l*s*exp(2*L*sqrt(s/D_l)) + sqrt(D_l*s)*(kH*sqrt(D_g*s) - sqrt(D_l*s))*(exp(2*L*sqrt(s/D_l)) + 1))*exp(2*l*sqrt(s/D_g))

                )
            )


        #0.5*sqrt(pi()/c)*(1-erf(s/(2*sqrt(pi()))))*exp(s**2/(4*c))
        #*1/(s+c)
        #1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))*

        kernel=np.zeros(len(x))
        for i in range(len(x)):
            kernel[i]=invertlaplace(fs,x[i],method='talbot')


        #disp = taylor_dispersion(x, D=c, L=l_cap)

        #kernel = signal.convolve(kernel, disp, mode='full')
        #kernel = kernel[:-(len(x)-1)]

        if matrix:
            kernel = np.tile(kernel,(len(kernel),1))
            i=1
            while i < len(x):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i=i+1

    if norm:
        area = np.trapz(kernel, x)
        return abs(kernel/area)
    else:
        return kernel

def gen_kernel4(x, D_l, L, A1, D_g, l, A2, kH, V_dot, fourier = False, matrix = False, norm=False, verbose = False):
    """
    Returns the kernel along x (time or frequency) with
        D: Diffusion coefficient
        L: Distance between WE and chip
        h: mass transport coefficient at liquid/gas interface
        c: pumping speed correction parameter
    In time domain the output can also be in form of the convolution matrix generated
    from the kernel.
    """
    if verbose:
        print('D: ' + str(D) + ' || L: ' + str(L) + ' || h: ' + str(h) + ' || c: ' + str(c))

    if fourier:
        print('Kernel in fourier space implemented soon')
        fs = lambda s: 1/(sqrt(s)*sinh(sqrt(s))+a*cosh(sqrt(s)))#*1/(2*c*sqrt(2*pi()))*exp(sqrt(s/c)/8)*erfc(s/(2*sqrt(2)*c))

    else:
        tdiff=x*D_g/(l**2)
        a = D_l*l/(D_g*L)
        b = l/L
        c = A2/A1
        Pe = V_dot*l/(A2*D_g)


        fs = (
            lambda s: -2.0*a*c*kH*sqrt(s/(a*b))*sqrt(Pe**2 + 4*s)*exp(1.0*Pe)/(
                sqrt(a*s/b)*(
                        (Pe - sqrt(Pe**2 + 4*s))
                        *(
                            a*c*sqrt(s/(a*b))*sinh(sqrt(s/(a*b)))
                            + 0.5*kH*(Pe - sqrt(Pe**2 + 4*s))*cosh(sqrt(s/(a*b)))
                            )
                        *exp(0.5*Pe - 0.5*sqrt(Pe**2 + 4*s))
                        - (Pe + sqrt(Pe**2 + 4*s))
                        *(
                            a*c*sqrt(s/(a*b))*sinh(sqrt(s/(a*b)))
                            + 0.5*kH*(Pe + sqrt(Pe**2 + 4*s))*cosh(sqrt(s/(a*b)))
                            )
                        *exp(0.5*Pe + 0.5*sqrt(Pe**2 + 4*s))
                    )

                )
            )

        kernel=np.zeros(len(x))
        for i in range(len(x)):
            kernel[i]=invertlaplace(fs,tdiff[i],method='talbot')


        #disp = taylor_dispersion(x, D=c, L=l_cap)

        #kernel = signal.convolve(kernel, disp, mode='full')
        #kernel = kernel[:-(len(x)-1)]

        if matrix:
            kernel = np.tile(kernel,(len(kernel),1))
            i=1
            while i < len(x):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i=i+1

    if norm:
        area = np.trapz(kernel, x)
        return abs(kernel/area)
    else:
        return kernel

def extract_kernel(dataset, tspan, mass, pot, t_bg = None):
    """
    Extracts a measured kernel of a defined mass at tspan from a dataset.
    """

    x_curr, y_curr = dataset.get_current(tspan=tspan)
    x_pot, y_pot = dataset.get_potential(tspan=tspan)
    x_sig, y_sig = dataset.get_signal(mass, tspan=tspan)

    if mass == 'M32':
        t0 = x_curr[np.argmax(y_pot > pot)] #time of current impulse
    elif mass == 'M2':
        t0 = x_curr[np.argmax(y_pot < pot)]  #time of current impulse
    else:
        print('mass not found')

    x_sig = x_sig - t0

    y_sig = y_sig[x_sig>0]
    x_sig = x_sig[x_sig>0]

    y_curr = y_curr[x_curr>t0]
    x_curr = x_curr[x_curr>t0]
    y_pot = y_pot[x_pot>t0]
    x_pot = x_pot[x_pot>t0]


    if t_bg is not None:
        _, y_bg = dataset.get_signal(mass, tspan=t_bg)
        y_sig = y_sig - np.mean(y_bg)

    kernel = Kernel(
        MS_data = np.array([x_sig, y_sig]),
        EC_data = np.array([x_curr, y_curr, x_pot, y_pot])
        )


    return kernel

def fit_kernel(dataset, tspans, mol, pot, initials = None, bounds = None, plot=False):
    """
    Fits measured kernel to analytical kernel extracting mass transport parameters.
        initials: initial guesses for the fitting parameters (see scipy_curvefit)
        bounds: bounds for fitting parames (see scipy_curvefit)
        plot: path + filename where output is saved as plot

    Returns:
        fitted L, h, c (see gen_kernel) and their standard deviation.

    #notes:
    Major speed-up probably through fitting in laplace domain instead of time domain
    """
    mol = Molecule(mol)
    if initials == None:
        #import initial guess from mol.py with 100um cell distance and 1.5e^15 molecules/s, high c initial means infinite pumping speed
        h = mol.kH * Chem.R * 1.5e15 / Chem.NA * 298 / (1e5*0.196e-4) #alot of hardcoded stuff -> not good xxxx
        initials = [100e-6, h, 6]

    if bounds == True:
        bounds=([50e-6,1e-10,1.5],[300e-6,1e-3,100])

    if any(isinstance(i, list) for i in tspans):
        x_kernel=[]
        y_kernel=[]
        for tspan in tspans:
            x,y = extract_kernel(dataset, tspan, mol.primary, pot, norm=True)

            x_kernel = np.concatenate((x_kernel, x))
            y_kernel = np.concatenate((y_kernel, y))
    else:
        x_kernel,y_kernel = extract_kernel(dataset, tspans, mol.primary, pot)
        y_kernel = y_kernel - y_kernel[0]
        y_kernel = y_kernel/np.max(y_kernel)
        x_kernel = x_kernel




    def kernel(x, L, h, c):
        return gen_kernel(x, mol.D, L, h, c, norm=True, verbose=True)

    pars, pcov = curve_fit(kernel, x_kernel, y_kernel, p0=initials, bounds=bounds)

    pstd = np.sqrt(np.diag(pcov))

    if plot:
        fig1 = plt.figure()
        axe1 = fig1.add_subplot(111)
        axe1.plot(x_kernel, y_kernel, linestyle = 'None', marker = '.')
        axe1.plot(x_kernel, gen_kernel(x_kernel, mol.D, pars[0], pars[1], pars[2], norm=True))
        fig1.savefig(mol.name + '_impulse' + str(plot) + '.png')

    return pars, pstd

def wiener_deconvolution(sig, kernel, lambd):
    "lambd is the inverse SNR"
    kernel = np.hstack((kernel, np.zeros(len(sig) - len(kernel)))) # zero pad the kernel to same length
    H = fft(kernel)
    deconvolved = np.real(ifft(fft(sig)*np.conj(H)/(H*np.conj(H) + lambd**2)))
    return deconvolved

def wiener_deconvolution_bg(sig, kernel, bg, freq):

    f, Pxx_noise = signal.welch(bg, freq)
    f, Pxx_sig = signal.welch(sig, freq)
    SNR = interp1d(f, Pxx_sig/(Pxx_sig+Pxx_noise))
    f_sig = fftfreq(len(sig), 1/freq)
    SNR_vec = SNR(abs(f_sig))
    kernel = np.hstack((kernel, np.zeros(len(sig) - len(kernel))))

    H = fft(kernel)
    SNR_vec = 1
    deconvolved = np.real(ifft(fft(sig)*np.conj(H)/(H*np.conj(H)+0.0001)*SNR_vec))

    return deconvolved

def get_rate(dataset, mass, tspan, t_bg, kernel_obj, lambd=0.1):

    x_sig, y_sig = dataset.get_signal(mass, tspan=tspan, t_bg=t_bg)

    time = np.arange(0,100,x_sig[1] - x_sig[0])
    time[0] = 1e-6
    kernel = kernel_obj.gen_kernel(time=time)

    y_sig = y_sig * sum(kernel)

    rate = wiener_deconvolution(y_sig, kernel, lambd)

    return x_sig, rate



class Deconvolution:
    """
    Class description.
    """
    def __init__ (
        self,
        D,
        L,
        h,
        c
    ):
        """
        All parameters in SI units.
        """
        self.length = L
        self.D = D
        self.h = h
        self.c = c


    @property
    def tau(self):
        return (self.length/self.h+self.length**2/(2*self.D)) #time constant

    def MS_to_EC(self, time, MS_signal, bg, kernel=None):

        freq = 1/abs(time[0]-time[1])
        k_time = np.concatenate(([1e-6], time[1:50]))

        #make that kernel can be passed directly
        if kernel is None:
            kernel = gen_kernel(k_time, self.D, self.length, self.h, self.c, norm=True)
        #kernel = kernel[kernel>0.1]
        #print(kernel)
        MS_signal=MS_signal * sum(kernel)

        while len(bg)<256:
            bg = np.concatenate((bg,bg))


        #MS_signal=gaussian_filter1d(MS_signal, 2)
        #MS_signal=signal.wiener(MS_signal,3)
        #MS_signal=gaussian_filter1d(MS_signal, 0.5)
        #MS_signal=signal.savgol_filter(MS_signal,9,3)
        #sig=interp1d(time-t_min, MS_signal)



        current = wiener_deconvolution_bg(MS_signal, kernel, bg, freq)


        #current = current-current[0]
        #zero_fill = np.zeros(2*t_length-1-len(current))
        #current = np.concatenate((current, zero_fill))
        #current = current*self.n_el*96485
        #current=gaussian_filter1d(current, 1)
        #current=signal.savgol_filter(current,3,1)
        #current=signal.wiener(current,mysize=10)
        #return current
        return current


    def EC_to_MS(self, time, EC_current):

        t_length = len(time)
        t_max = time[-1]
        t_min = time[0]


        time_lin, step = np.linspace(0,t_max-t_min, t_length, retstep=True) #converts the input time to a linear array to ensure evenly spaced time values


        curr=interp1d(time-t_min, abs(EC_current/(self.n_el*96485)))
        impulse=self.trans_num(time_lin)
        conc_interphase=signal.convolve(curr(time_lin), impulse, mode='full') / sum(trans_num(time_lin)) * self.length/self.D

        flux_MS = conc_interphase*self.h
        p_i = flux_MS/self.capillary_flow
        p_i = interp1d(np.linspace(0,2*(t_max-t_min)-step, 2*t_length-1), p_i)

        return p_i(time-t_min)

class Kernel:
    """
    Class description.
    """
    def __init__ (
        self,
        parameters = {},
        MS_data = None,
        EC_data = None,
    ):
        """
        All parameters in SI units.
        """
        if MS_data is not None and parameters:
            raise Exception(
                'Kernel can only be initialized with data OR parameters, not both'
                )
        if EC_data is not None and MS_data is not None:
            print('Generating kernel from measured data')
            self.type = 'measured'
        elif parameters:
            print('Generating kernel from parameters')
            self.type = 'functional'
        else:
            print('Generating blank kernel')
            self.type = None

        self.params = parameters
        self.MS_data = MS_data
        self.EC_data = EC_data #x_curr, y_curr, x_pot, y_pot

    @property
    def sig_area(self):

        delta_sig = self.MS_data[1] - self.MS_data[1][-1]
        sig_area = np.trapz(delta_sig, self.MS_data[0])

        return sig_area


    @property
    def charge(self):

        y_curr = self.EC_data[1]

        mask = np.isclose(
            y_curr,
            y_curr[0],
            rtol=1e-1
            )

        Q = np.trapz(y_curr[mask], self.EC_data[0][mask])

        return Q

    def plot(self,
        time=None,
        ax = None,
        norm=True,
        **kwargs
        ):

        if ax is None:
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)

        if self.type is 'functional':
            ax.plot(
                time, self.gen_kernel(time=time, norm=norm),
                **kwargs,
                )

        elif self.type is 'measured':
            ax.plot(
            self.MS_data[0], self.gen_kernel(norm=norm),
            **kwargs,
            )

        else:
            raise Exception('Nothing to plot with blank kernel')


        return ax

    def gen_kernel(self,
        time=None,
        norm=True,
        matrix=False
        ):

        if self.type is 'functional':

            if time is None:
                raise Exception('Functional kernel generation requires time kwarg')

            D = self.params['D']
            L = self.params['L']
            V = self.params['V']
            V_dot = self.params['V_dot']
            kH = self.params['kH']
            c = self.params['c']

            tdiff = time * D / (L**2)
            fs = lambda s: 1/(
                sqrt(s) * sinh(sqrt(s))
                + (V * kH / 0.196e-4 / L)
                * (s + V_dot / V * L**2 / D)
                * cosh(sqrt(s))
                )*1/(s+1/c*L**2/D)

            kernel = np.zeros(len(time))
            for i in range(len(time)):
                kernel[i] = invertlaplace(fs, tdiff[i], method='talbot')

            #dx = time[1] - time[0]
            #gx = np.arange(-7*c, 7*c, dx)
            #gaussian = np.exp(-(gx/c)**2/2)
            #kernel = np.convolve(kernel, gaussian, mode="full") / np.sum(gaussian)
            #kernel = kernel[0:len(time)]

        elif self.type is 'measured':
            kernel = self.MS_data[1]
            time = self.MS_data[0]

        if norm:
            area = np.trapz(kernel, time)
            kernel = kernel/area

        if matrix:
            kernel = np.tile(kernel, (len(kernel), 1))
            i = 1
            while i < len(x):
                kernel[i] = np.concatenate((kernel[0][i:], kernel[0][:i]))
                i=i+1


        return kernel

    def fit(self, mol = 'H2', initials=None, bounds=None):
        #under construction xxxx
        mol = Molecule(mol)

        def kernel(x, L, V, V_dot):
            params = {
                'D' : mol.D,
                'L' : L,
                'V_dot' : V_dot,
                'V' : V,
                'kH' : mol.kH,
                }

            fit_kernel = Kernel(params=initials)

            return fit_kernel.gen_kernel()

        pars, pcov = curve_fit(kernel, x_kernel, y_kernel, p0=initials, bounds=bounds)

        pstd = np.sqrt(np.diag(pcov))

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:28:25 2020

@author: scott
"""

import os, re
import numpy as np
from matplotlib import pyplot as plt

data_directory = os.path.dirname(os.path.realpath(__file__)) + os.sep + "data"




def gauss(x, center, sigma, height):
    y = height * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))
    return y


class Peak:
    """
    Copied from https://github.com/ScottSoren/EC_Xray/blob/master/XRD.py
    on April 1, 2020.
    """

    def __init__(self, x, y, xspan=None, name=None, color="k"):
        self.name = name
        self.xspan = xspan
        self.color = color
        if xspan is not None:
            mask = np.logical_and(xspan[0] < x, x < xspan[-1])
            x, y = x[mask], y[mask]
        self.x = x
        self.y = y
        self.background = np.zeros(x.shape)
        self.bg = False

    def set_background(self, background):
        if not len(background) == len(self.x):
            print("background wrong length.")
            raise ValueError
        self.background = background
        self.bg = True

    def get_background(self, *args, **kwargs):
        x, y = self.x, self.y
        self.background = get_peak_background(x, y, **kwargs)
        self.bg = True
        return self.background

    def get_integral(self, *args, **kwargs):
        if "mode" in kwargs and kwargs["mode"] in ["gauss", "fit"]:
            if "ax" in kwargs:
                self.fit_gauss(ax=kwargs["ax"])
            return self.integral_f
        x, y, background = self.x, self.y, self.background
        integral = np.trapz(y - background, x)
        self.integral = integral
        if "ax" in kwargs:
            ax = kwargs["ax"]
            if ax == "new":
                fig, ax = plt.subplots()
            if ax is not None:
                ax.plot(x, y, "k.")
                ax.plot(x, background, "b--")
                ax.fill_between(x, background, y, where=y > background, color="g")
        return integral

    def fit_gauss(self, center=None, sigma=None, ax=None):
        x, y, background = self.x, self.y, self.background
        y = y - background

        guess_c = (x[-1] + x[0]) / 2
        guess_s = (x[-1] - x[0]) / 2
        guess_h = max(y)

        if center is not None and sigma is not None:

            def gauss_i(x, height):
                return gauss(x, center=center, sigma=sigma, height=height)

            guess = guess_h
            popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
            height = popt[0]
        elif center is not None:

            def gauss_i(x, sigma, height):
                return gauss(x, center=center, sigma=sigma, height=height)

            guess = [guess_s, guess_h]
            popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
            sigma, height = popt[0], popt[1]
        elif sigma is not None:

            def gauss_i(x, center, height):
                return gauss(x, center=center, sigma=sigma, height=height)

            guess = [guess_c, guess_h]
            popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
            center, height = popt[0], popt[1]
        else:

            def gauss_i(x, center, sigma, height):
                return gauss(x, center=center, sigma=sigma, height=height)

            guess = [guess_c, guess_s, guess_h]
            try:
                popt, pcov = curve_fit(gauss_i, x, y, p0=guess)
                center, sigma, height = popt[0], popt[1], popt[2]
            except RuntimeError:
                center, sigma, height = guess
        sigma = abs(sigma)
        # print(f'center={center}, sigma={sigma}, height={height}') # debugging
        fit = gauss(x, center, sigma, height)
        integral_f = np.sqrt(2 * np.pi) * height * sigma
        self.center, self.sigma, self.height = center, sigma, height
        self.fit, self.integral_f = fit, integral_f

        if ax is not None:
            if ax == "new":
                fig, ax = plt.subplots()
            ax.plot(x, background, "b--")
            ax.plot(x, y + background, "k.")
            ax.plot(x, fit + background, "r--")

        return center, sigma, height

    def get_fwhm(self):
        fwhm = 2 * np.sqrt(2 * np.log(2)) * self.sigma
        return fwhm

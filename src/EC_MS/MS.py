# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:28:25 2020

Peak class copied from EC_Xray commit 13b3dbf

@author: scott
"""

import os, re
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

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
        self.y_bg = np.zeros(x.shape)
        self.bg = False

    def get_background(self, y_bg=None, bg_mode="linear", endpoints=2):
        """
        This can be much simpler than the ridiculous background subtraction
        stuff I had going in EC_Xray because mass spec peaks always come
        at integer m/z and so it's easy to surround the peak
        """
        x, y = self.x, self.y
        if bg_mode == "linear":
            # draw a background line through the endpoints
            x_start, x_finish = np.mean(x[:endpoints]), np.mean(x[-endpoints:])
            y_start, y_finish = np.mean(y[:endpoints]), np.mean(y[-endpoints:])

            y_bg = y_start + (x - x_start) / (x_finish - x_start) * (y_finish - y_start)

        return y_bg

    def set_background(self, bg=True, y_bg=None, bg_mode="linear", endpoints=2):
        self.reset()
        if bg:
            y_bg = self.get_background(y_bg=y_bg, bg_mode=bg_mode, endpoints=endpoints)
            self.y_bg = y_bg
        else:
            self.y_bg = np.zeros(self.x.shape)

    def subtract_background(self, y_bg=None, bg_mode="linear", endpoints=2):
        self.reset()  # to avoid losing the ability to restore the original
        # by subtracting a new background from background-subtracted data
        y_bg = self.get_bacground(y_bg=y_bg, bg_mode=bg_mode, endpoints=endpoints)

        self.y_bg = y_bg
        self.bg = True
        self.y = self.y - y_bg

    def reset(self):
        """
        so far only implemented for masses.
        """
        if self.bg:
            y_bg = self.y_bg
            self.y = self.y + y_bg
            self.bg = False

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

    def fit_gauss(
        self, center=None, sigma=None, ax=None, y_bg=None, bg_mode=None, endpoints=2
    ):
        if y_bg or bg_mode:
            self.set_background(y_bg=y_bg, bg_mode=bg_mode, endpoints=endpoints)

        x, y, y_bg = self.x, self.y, self.y_bg
        y = y - y_bg

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
            ax.plot(x, y_bg, "b--")
            ax.plot(x, y + y_bg, "k.")
            ax.plot(x, fit + y_bg, "r--")

        return center, sigma, height

    def get_fwhm(self):
        fwhm = 2 * np.sqrt(2 * np.log(2)) * self.sigma
        return fwhm

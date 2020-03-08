# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:55:48 2017

@author: murems
"""

import os, codecs, io
import numpy as np

U_lines = ["Ei (V)", "E1 (V)", "E2 (V)", "Ef (V)", "EL (V)"]
check_lines = ["vs."]


class SettingsFile:
    def __init__(
        self, name, file_name=None, RE_vs_RHE=None, settings={}, verbose=True, **kwargs
    ):
        settings.update(kwargs)
        self.settings = settings
        self.RE_vs_RHE = RE_vs_RHE
        if file_name is None and len(settings) == 0:
            file_name = name
        self.name = name
        self.verbose = verbose
        try:
            with open(file_name) as f:
                self.lines = f.readlines()
        except UnicodeDecodeError:
            print(
                "why on earth don't you idiots just use utf-8 for everything? Trying iso8859_15."
            )
            with codecs.open(file_name, encoding="iso8859_15") as f:
                self.lines = f.readlines()

    def write(self, file_name=None, encoding="iso8859_15"):
        if file_name is None:
            file_name = self.name + ".mps"
        with io.open(file_name, "w", encoding=encoding) as f:
            f.writelines(self.lines)

    def change_potentials(self, RE_vs_RHE=None, diff=None, verbose=None):
        if verbose is None:
            verbose = self.verbose
        if verbose:
            print("\n\nfunction 'change_potentials' at your service!\n")

        if RE_vs_RHE is not None:
            diff = -(RE_vs_RHE - self.RE_vs_RHE)
            self.RE_vs_RHE = RE_vs_RHE
        newlines = []
        N_tech = 0
        newtech = False
        N_lines = len(self.lines)
        for n in range(
            N_lines
        ):  # len(self.lines changes each step when iterating over inumerate. Wierd.)
            line = self.lines[n]
            if line[0:9] == "Technique":
                N_tech += 1
                if verbose:
                    print("Technique " + str(N_tech))
                newtech = True
                newline = line
            elif newtech:  # newtech stays true all through techniques like OCV
                if verbose:
                    print(line)
                newline = line
                if len(line) > 2 and line[0:2] == "Ns":
                    N_step = (len(line)) / 20 - 1
                    # Hmmm, the number of characters per line is different based on OS.
                    # Should be len(line)-2 in ubuntu. But I'll just round down later anyway.
                    if verbose:
                        print(str(N_step) + " steps.")
                    newtech = False  # now we're into the numerical settings!
                # but there is no Ns if the technique only has one step!
                elif (
                    n < N_lines
                    and len(self.lines[n + 1]) > 5
                    and self.lines[n + 1][0:6] in U_lines
                ):
                    N_step = 1
                    if verbose:
                        print("1 step")
                    newtech = False
            elif N_tech == 0:  # Header ends up here.
                newline = line
            else:
                parts = [
                    line[20 * i : 20 * (i + 1)] for i in range(int(N_step + 1))
                ]  # +1 for line title
                #                print(parts)
                if len(parts) == 0 or n == N_lines:
                    if verbose:
                        print("skipping line " + str(n) + " out of " + str(N_lines))
                    newline = line
                else:
                    if verbose:
                        print("row: " + parts[0])
                    if parts[0][0:6] in U_lines:
                        check_line = self.lines[n + 1]
                        check_parts = [
                            check_line[20 * i : 20 * (i + 1)]
                            for i in range(int(N_step + 1))
                        ]  # +1 for line title
                        newline = ""
                        for part, check_part in zip(parts, check_parts):
                            print(part + " " + check_part)
                            if check_part[0:3] == "Ref":
                                U = eval(part)
                                U_new = np.round(U + diff, 3)
                                if verbose:
                                    print(
                                        "replacing "
                                        + str(U)
                                        + " V vs Ref with "
                                        + str(U_new)
                                        + " V vs Ref."
                                    )
                                newpart = str(U_new).ljust(20)
                                newline += newpart
                            else:  # line title also ends up here.
                                newline += part
                        newline += "\n"
                    else:
                        newline = line
            newlines += [newline]
        self.lines = newlines


if __name__ == "__main__":
    data_dir = os.path.expanduser(
        "~/Dropbox/Sniffer_Experiments/06_HER/Data/17G10_Nik3/"
    )
    os.chdir(data_dir)
    test = SettingsFile(
        name="test", file_name="closing_secondary.mps", RE_vs_RHE=1.298, verbose=True
    )
    test.change_potentials(RE_vs_RHE=0.5900)
    test.write(file_name="RE_vs_RHE_is_0p59V.mps")

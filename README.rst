=====================================================================
``EC_MS``: Electrochemistry plus Mass Spectrometry
=====================================================================

``EC_MS`` provides a powerful **object-oriented** interface to electrochemistry data, mass spectrometry data,
and especially the combination of these two types of datasets produced by electrochemistry - mass spectrometry (EC-MS) techniques such as
differential electrochemical mass spectrometry (DEMS) and chip-based EC-MS.

``EC_MS`` has grown in concert with the chip EC-MS technology sold by `Spectro Inlets <https://spectroinlets.com>`_, but supports analysis of data from other hardware.

The primary object-oriented interface for this is the ``Dataset`` class. For example:

.. -code-begin-
.. code-block:: pycon

   >>> from EC_MS import Dataset
   >>> MS_dataset = Dataset('MS_data.txt', data_type='MS')
   >>> EC_dataset = Dataset('EC_data.mpt', data_type='EC')
   >>> dataset = MS_dataset + EC_dataset # calls the function EC_MS.synchronize()
   >>> dataset.plot_experiment() # EC data in lower panel, MS data in upper panel

In this example, the MS and EC datasets are combined by lining up all of the time variables based on timestamps read in the headers of the files.

It is easy to manipulate the datasets based on the electrochemistry program

.. -code-begin-
.. code-block:: pycon

   >>> from EC_MS import CyclicVoltammagram
   >>> cv = CyclicVoltammagram(Dataset)
   >>> cv.normalize(RE_vs_RHE=0.715)
   >>> cv.redefine_cycle(V=0.45, redox=1) # defines when the cycle counter increases
   >>> cycle_1 = cv[1] # selects one cycle
   >>> cycle_1.plot(masses=['M2', 'M44']) # electrochemical potential on the x-axis

And that's just a small teaser. Additional functionality includes:

- object-oriented interface to mass spectra with the Spectrum and Spectra classes

- Calibration functions and classes for quantitative data analysis and plotting

- Thermochemistry and Electrolyte subpackages for calculating standard potentials and chemical equilibrium

- Mass-transport modelling of products and reactants in the working volume between the electrode and the vacuum inlet

- ohmic drop correction and automated quantitative comparisons of cyclic voltammagrams


Full documentation is pending!


Installation
============

EC_MS is pip-installable! Just type in your terminal or Anaconda prompt:

.. -code-begin-
.. code-block:: bash

   $ pip install EC_MS

The in-development version is available on `github <https://github.com/ScottSoren/EC_MS/>`_.

**EC_MS** requires **numpy**, **scipy**, and **matplotlib**. I recommend using *Anaconda* python, and writing and running your scripts with *spyder*. This has proven the easiest to set up on all operating systems I've tried.


Supported Data Types
====================

**Mass Spectrometry**

- .tsv files from Spectro Inlets' Zilien (data_type="SI")

- .dat files (both Bin.dat and Scan.dat) from Pfeiffer Vacuum's PVMassSpec (data_type="PVMS")

- .txt files from `cinfdata <https://github.com/CINF/cinfdata>`_. (data_type="MS")

- .txt files from Stanford Reasearch Systsms' Residual Gas Analyzer (data_type="RGA")


**Electrochemistry**

- .tsv files from Spectro Inlets' Zilien (data_type="SI")

- .mpt files from BioLogic's EC-Lab (data_type="EC")

- .txt files from CH Instruments software (data_type="CHI")

Full documentation is pending!


If you would like support for another file type, write to me.


References
==========

This python package was first described in:

Daniel B. Trimarco and Soren B. Scott, et al. **Enabling real-time detection of electrochemical desorption phenomena with sub-monolayer sensitivity**. `Electrochimica Acta, 2018 <https://doi.org/10.1016/j.electacta.2018.02.060>`_.

Its functionality is demonstrated, a bit more up-to-date, in the figures and footnotes of:

Soren B. Scott. **Isotope-Labeling Studies in Electrocatalysis for Renewable Energy Conversion and the Net CO2 Impact of this PhD Project.** `PhD Thesis, 2019. <https://orbit.dtu.dk/en/publications/isotope-labeling-studies-in-electrocatalysis-for-renewable-energy>`_.

Other articles with figures and data analysis by **EC_MS** include:

- Anna Winiwarter and Luca Silvioli, et al. Towards an Atomistic Understanding of Electrocatalytic Partial Hydrocarbon Oxidation: Propene on Palladium. `Energy and Environmental Science, 2019 <https://doi.org/10.1039/C8EE03426E>`_.

- Claudie Roy, Bela Sebok, Soren B. Scott, et al.  Impact of nanoparticle size and lattice oxygen on water oxidation on NiFeOxHy. `Nature Catalysis, 2018 <https://doi.org/10.1038/s41929-018-0162-x>`_.



Project Information
===================

This is a pre-alpha version, so it is buggy. Please log issues on `github <https://github.com/ScottSoren/EC_MS/>`_ to help me improve it.

``EC_MS`` is fully free and open-source.

If you have questions or if you'd like to contribute, please write to me.
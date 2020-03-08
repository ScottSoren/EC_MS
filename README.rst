=====================================================================
``EC_MS``: Electrochemistry and Mass Spectrometry
=====================================================================

``EC_MS`` provides a powerful **object-oriented** interface to electrochemistry data, mass spectrometry data, and especially the combination of these two types of datasets produced by electrochemistry - mass spectrometry techniques such as differential electrochemical mass spectrometry (DEMS) chip-based electrochemistry mass spectrometry (chip EC-MS). 

The primary object-oriented interface for this is the ``Dataset`` class. For example:

.. -code-begin-
.. code-block:: pycon

   >>> from EC_MS import Dataset
   >>> MS_dataset = Dataset('MS_data.txt', data_type='MS')
   >>> EC_dataset = Dataset('EC_data.mpt', data_type='EC')
   >>> dataset = MS_dataset + EC_dataset
   >>> dataset.plot_experiment()

In this example, the MS and EC datasets are combined by lining up all of the time variables based on timestamps read in the headers of the files. 

It is easy to manipulate the datasets based on the electrochemistry program

.. -code-begin-
.. code-block:: pycon

   >>> from EC_MS import CyclicVoltammagram
   >>> cv = CyclicVoltammagram(Dataset)
   >>> cv.normalize(RE_vs_RHE=0.715)
   >>> cv.redefine_cycle(V=0.45, redox=1) # defines when the cycle counter increases
   >>> cycle_1 = cv[1]
   >>> cycle_1.plot(masses=['M2', 'M44'])

And that's just a small teaser. Additional functionality includes: 
- Calibration functions and classes for quantitative data analysis and plotting
- Thermochemistry and Electrolyte subpackages for calculating standard potentials and chemical equilibrium
- Mass-transport modelling of products and reactants in the working volume between the electrode and the vacuum inlet
- ohmic drop correction and automated quantitative comparisons of cyclic voltammagrams

References
==========

This python package was first described in:

Daniel B. Trimarco and Soren B. Scott, et al. **Enabling real-time detection of electrochemical desorption phenomena with sub-monolayer sensitivity**. `Electrochimica Acta, 2018 <https://doi.org/10.1016/j.electacta.2018.02.060>`_.

The theory behind its quantification tools are described more thoroughly and correctly in:

Soren B. Scott. **Isotope-Labeling Studies in Electrocatalysis for Renewable Energy Conversion and the Net CO2 Impact of this PhD Project.** `PhD Thesis, 2019. <https://orbit.dtu.dk/en/publications/isotope-labeling-studies-in-electrocatalysis-for-renewable-energy>`_.

Project Information
===================

If you have questions or if you'd like to contribute, please log issues on `github <https://github.com/ScottSoren/EC_MS/>`_ or write to me.
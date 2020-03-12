.. EC_MS documentation master file, created by
   sphinx-quickstart on Thu Mar 12 20:35:48 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EC_MS's documentation!
=================================

EC_MS is the python module that will really "super-charge" your EC-MS
data treatment and plotting.

.. code-block:: pycon

   >>> from EC_MS import Dataset
   >>> MS_dataset = Dataset('MS_data.txt', data_type='MS')
   >>> EC_dataset = Dataset('EC_data.mpt', data_type='EC')
   >>> dataset = MS_dataset + EC_dataset  # calls the function EC_MS.synchronize()
   >>> dataset.plot_experiment()          # EC data in lower panel, MS data in upper panel

To get started with your EC-MS data treatment in a hurry, read the
:ref:`getting started <getting_started>` page, with :ref:`installation
instructions <installation>` and a small :ref:`tutorial <tutorial>`.

To search the complete code documentation head over to the :ref:`API
docs <module_reference>`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api/EC_MS.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

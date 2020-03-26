.. _getting_started:

***************
Getting started
***************

This section will help you to quickly get started with *EC_MS*.

.. _installation:

Installation
============

*EC_MS* can be installed either with :ref:`pip <installation_pip>`
(recommended) or :ref:`manually <installation_manually>`.

.. _installation_pip:

From PyPi with pip
------------------

The easiest way to install *EC_MS*, is to install it from `PyPi
<https://pypi.python.org/pypi>`_ with the program `pip
<https://pip.pypa.io/en/stable/>`_. This can be done with the command:

.. code-block:: sh

   pip install EC_MS

This will automatically take care of installing any dependencies you
need.

.. _installation_manually:

Manual installation from .tar.gz file
-------------------------------------

*EC_MS* can also be installed manually from the .tar.gz file. First,
find `the latest version of EC_MS on PyPi
<https://pypi.org/project/EC-MS/#files>`_ and download the
newest ``.tar.gz`` file. After that, extract the content and move into
the extracted folder. As an example, for *EC_MS* 0.5.1 and on a Unix
type system, this can be done with the following commands:

.. code-block:: sh

   wget https://files.pythonhosted.org/packages/d3/02/bb8cc9dcc70def1339843c6ef8d3dcc02cf5e22cbb9925097c9ab2b36d5f/EC_MS-0.5.1.tar.gz
   tar zxvf EC_MS-0.5.1.tar.gz
   cd cd EC_MS-0.5.1/

Have a look inside the ``requirements.txt`` file. You will need to install the
dependencies listed in that file yourself. See the documentation for the
individual dependencies for installation instructions.

After the requirements are in place, the package can be install with the
command:

.. code-block:: sh

   python setup.py install

After installation check
------------------------

After installation, open a Python interpreter and check that :mod:`EC_MS` can
be imported:

.. code-block:: python

   >>> import EC_MS

.. _tutorial:

Tutorial
========

TODO

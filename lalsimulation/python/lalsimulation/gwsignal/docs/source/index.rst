.. GWSignal documentation master file, created by
   sphinx-quickstart on Wed Jul  6 10:06:13 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GWSignal documentation!
==================================


GWSignal
----------------
:code:`GWSignal` is a Python library that provides a interface between python waveforms and `LALSuite <https://lscsoft.docs.ligo.org/lalsuite/>`_ and vice-versa. This software is packaged with  
`LALSimulation <https://lscsoft.docs.ligo.org/lalsuite/>`_ in `LALSuite <https://lscsoft.docs.ligo.org/lalsuite/>`_ releases. Check out the :doc:`installation` section for further information, including how to 
:ref:`install <install>` the project.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   overview
   usage
   implementing_your_model
   parameters_inputs
   gw_parameters
   examples


.. toctree::
   :maxdepth: 1
   :caption: Examples:

   examples/example_usage.ipynb
   examples/testing_conditioning_routines.ipynb   
   examples/error_handling.ipynb   

	

API
---
.. autosummary::
   :toctree: api
   :template: custom-module-template.rst
   :caption: API:
   :recursive:

   gwsignal.core
   gwsignal.models

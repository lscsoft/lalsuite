Installation
============
Quick guide for installing :code:`GWSignal` to your system.

.. _install:

Install GWSignal with Pip
---------------------------
As :code:`GWSignal` is packaged with LALSuite, it can be pip installed from the `lalsuite <https://pypi.org/project/lalsuite/>`_ pip package.

.. code-block:: console

   (.home) $ pip install lalsuite


Install GWSignal from Source
------------------------------

Clone the LALSuite `repository <https://git.ligo.org/lscsoft/lalsuite>`_ and install the software :

.. code-block:: console

    (.home) $ git clone git@git.ligo.org:lscsoft/lalsuite.git
    (.home) $ conda env create -f common/conda/environment.yml
    (.home) $ conda activate lalsuite-dev
    (.lalsuite-dev) $ ./00boot
    (.lalsuite-dev) $ ./configure --prefix=${CONDA_PREFIX}
    (.lalsuite-dev) $ make
    (.lalsuite-dev) $ make install

See the LALSuite `repository <https://git.ligo.org/lscsoft/lalsuite>`_ for further details / options for the :code:`./configure` step.

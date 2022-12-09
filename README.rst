#################
Simple IBM Solver
#################

THIS IS AN ON-GOING PROJECT.

|License|_ |CI|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleIBMSolver
.. _License: https://opensource.org/licenses/MIT

.. |CI| image:: https://github.com/NaokiHori/SimpleIBMSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleIBMSolver/actions/workflows/ci.yml

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleIBMSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleIBMSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleIBMSolver/blob/main/docs/source/snapshot2d.png
   :width: 100%

.. image:: https://github.com/NaokiHori/SimpleIBMSolver/blob/main/docs/source/snapshot3d.png
   :width: 50%

********
Overview
********

This library numerically solves the motion of moving rigid particles in two- and three-dimensional Cartesian domains using finite-difference and immersed boundary methods.
This is built on top of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_

Docker image is available.

.. code-block:: console

   $ mkdir /path/to/your/working/directory
   $ cd    /path/to/your/working/directory
   $ git clone https://github.com/NaokiHori/SimpleIBMSolver
   $ cd SimpleIBMSolver
   $ docker build -t simplenavierstokessolver:latest .

***********
Quick start
***********

Please check the README of `SimpleNavierStokesSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.

****************
Acknowledgements
****************

The immersed boundary method is based on `a publication <https://www.sciencedirect.com/science/article/pii/S0045793021003716>`_ with quite a few modifications.


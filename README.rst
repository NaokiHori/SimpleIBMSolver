#################
Simple IBM Solver
#################

|License|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleIBMSolver
.. _License: https://opensource.org/licenses/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleIBMSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleIBMSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleIBMSolver/blob/main/docs/source/thumbnail.gif
   :width: 100%

.. image:: https://github.com/NaokiHori/SimpleIBMSolver/blob/main/docs/source/snapshot3d.png
   :width: 100%

********
Overview
********

This library numerically solves the motion of moving rigid particles in two- and three-dimensional Cartesian domains using the finite-difference and the immersed boundary methods.

**********
Dependency
**********

This is built on top of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNSSolver>`_.
Please visit it and check the dependency there.

***********
Quick start
***********

#. Prepare workplace

   .. code-block:: console

      mkdir -p /path/to/your/directory
      cd       /path/to/your/directory

#. Get source

   .. code-block:: console

      git clone --recurse-submodules https://github.com/NaokiHori/SimpleIBMSolver
      cd SimpleIBMSolver

#. Set initial condition

   Here ``Python3`` is used to initialise the flow fields conveniently.
   One can give ``NPY`` files in different way under ``initial_condition/output/``.

   .. code-block:: console

      cd initial_condition
      make output
      bash exec.sh
      cd ..

#. Build solver

   .. code-block:: console

      make output
      make all

#. Run

.. code-block:: console

   bash exec.sh

****
Note
****

The immersed boundary method and the collision model is based on `a publication <https://www.sciencedirect.com/science/article/pii/S0045793021003716>`_ with some modifications.


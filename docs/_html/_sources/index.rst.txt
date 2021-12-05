.. APySPAM documentation master file, created by
   sphinx-quickstart on Thu Nov 25 19:59:38 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to APySPAM's documentation!
===================================

``APySPAM`` is a Python package for running interacting galaxy simulations. It takes in 16 different parameters, and returns the morphology, star formation rate distribution and flux distribution of the resultant systems.

What it can do
---------------

The Advanced Python Stellar Animation Module (APySPAM) is a numerical galaxy interaction modeller based on the earlier work of John Wallin and Anthony Holincheck in the Java Stellar Particle Animation Module (JSPAM). From a selection of parameters that the user 
can tweak, this algorithm will:

   - Simulate a close flyby of two galaxies with the given parameters;
   - Calculate the expected change in flux distribution based on the interaction;
   - Accurately recreate the morphology of the system;
   - Save a flux distribution map of the resultant system;
   - Save a .txt file of the particle positions, velocities, colour-band fluxes and final star formation rates;
   - Save a .npy file of the particle positions, velocities, colour-band fluxes and final star formation rates.

For a full description of the astrophysics of the algorithm, as well as some further use cases please see the output paper that came with this documentation (ArXiv:$INSERT_ARCHIVE_NO). The aim of this documentation is to give a very high level overview of each of the 
different scripts the user might have to interct with, some example results as well as some basic and advanced tutorials in its use. I hope that the user finds this algorithm helpful with the work they are pursuing! Any issues or errors that you encounter in the
algorithm, please contact the author this documentation David O'Ryan at d.oryan@lancaster.ac.uk.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation/init
   quickstart/init
   script_documentation/init
   advanced_examples/init
   results/init
   example_use_case/emcee_example/init



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

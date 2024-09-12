# ts-mcmc
This repo is a simple prototype for using tskit to store phylogenetic MCMC output.
It currently features a simple MCMC algorithm targeting the Kingman coalescent without mutation.

The project contains the following:

- mcmc.py: The specification of the MCMC loop. Operates on a tree specified in `tree.py`, with successive iterates compactly stored into a tree sequence by a recorder class specified in `recorder.py`.
- recorder.py: The class for efficiently copying MCMC iterates (stored as tree objects specified in `tree.py`) into a tree sequence.
- run.py: A short script for running one-off simulations.
- tree.py: The data structure for storing and updating an individual tree in the MCMC loop. Also contains the specifications of the Kingman coalescent target distribution, and the MCMC proposal distribution.
- verification.py: A statistical comparison of TMRCAs, total branch lengths, and youngest clade probabilities against analytical means. Set to $10^5$ MCMC iterations by default, with a runtime of around 6 minutes. Runtime is roughly linear in the number of iterations.

In addition, the `C` folder provides implementations of all of the above in C.
In the `C` folder, call `make run` to compile the MCMC routine and `make verification` to compile the statistical comparison of simulation output to analytical values.
To run verification after compilation, navigate to the `C` folder, call `./verification > verification.txt` to produce and store simulation output, and then run `verification.R` to produce plots.
The C code is much faster than Python, and a call to `./verification` should complete in around 15 seconds for the default $10^5$ iterations.

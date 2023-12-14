#RK_tf2

The present code is intended to be used with the Heat and Mass Transfer
Technological Centre (CTTC) in-house code TermoFluids Algebraic (TFA), in which
a general implementation for the numerical time integration of the
Navier-Stokes equations is performed by means of the set Runge-Kutta schemes,
following the paper "Accuracy analysis of explicit Runge-Kutta methods applied
to the incompressible Navier-Stokes equations" by B. Sanderse, and B. Koren.
For a description of the method, please refer to this paper.

## Authors

The main structure of the code has been developed by Josep Plana-Riu, of the
Heat and Mass Transfer Technological Centre (CTTC), C/ Colom 11, 08222
Terrassa, Spain; which includes the implementation and validation of the
Runge-Kutta schemes in the main structure of the in-house code, mainly
implemented by Guillem Colomer, of CTTC.

## Repository structure

The repository is organized as follows. The 'tests' directory contains the MMS
simulations ran to verify the implementation of the Runge-Kutta as well as some
tests performed to check the conditions for monotonicity and phase-preserving
of the solutions. The 'mods' directory contains all the required functions to
integrate the Navier-Stokes equations with Runge-Kutta with a self-adaptive
time-step standpoint (in the 'sat' directory); the energy equation block is
still not tested. The directories 'channel_flow' and 'rayleigh_benard' contain
the configuration files as well as the specific functions to set up the case in
the directory 'specific_mods', as well as some additional functions to
represent the results. Eventually, the directory 'parallel_in_time' contains
the set of simulation directories (of the same kind of the 'channel_flow' and
'rayleigh_benard' cases) set up to be run with multiple simulations at a time.

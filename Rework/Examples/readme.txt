The attached scripts are version 0 of the Cambridge University Rapid Turbomachine Inverse System (CURTIS).

The method is outlined in Clark 2019 - "A step towards an intelligent aerodynamic blade design system" - ASME Turbo expo

The general idea is that by doing lots of inverse calcs at once their cost can be reduced compared to doing them individually at the time of need. This was exploited to generate a large (~2000) database of blade geometries that fullfill given loading styles and duties. The attached scripts then perform interpolation using RBF networks between these entries to enable real time blade section generation and loss prediction.

The current configuration is operational on python 2.6.6 and uses the following libraries:
-numpy
-scipy
-matplotlib
-pickle
-os
as well as a few more minor libraries.

The key components in the folder are:
CURTIS_V0.py - main script including GUI and calling to other functions  --- RUN THIS.
CURTIS_FUNCTIONS_V0.py - contains many of the subroutines required
/Database - contains the raw database for regeneration of the RBF networks if required
/Pickles - contains pickled versions of the trained networks

There are also some mises scripts that have been modified to suppress plotting routines and some shell scripts used to call them.
It is likeley that slight changes are required to get mises to communicate to the python on different systems or different mises versions.
However, it is hoped that the interpolation setup should work fairly standalone without too many modifications.

In CURTIS_V0.py @ line 686: load = 1. This causes the RBF networks to be fitted - their stored versions are quite large (~100MB).
After this first run set load = 0 and they will be loaded much faster from the stored directory.

Any queries can be directed to me at cjc95@cam.ac.uk

All best,
Chris


#####################################################################################

Version 0 notes:

- Contains the fundamental parameterisation and interpolation routines as well as the proof of concept database aimed at lowspeed HP turbine profiles.

- initial database used for paper:
-- varies inlet angle (-50 -> -20), exit angle (65 -> 72.5), exit Mach (0.4 -> 0.7)
-- varies Mpeak (1.2 -> 1.4), Lpeak (0.5 -> 0.6), Msss (1.2->1.8) and Mps(0.8->1.0)

-Distributed database does NOT vary:
--Gamma - gas composition (only matters at high Mach)
--AVDR  - axial velocity density ratio
--GTE   - Extra style variable
--Re    - Reynolds number is fixed at 600,000

- it does not include "Mode Switching" to prescribe pitch-to-chord rather than one of the aerodynamic variables

- it does not include Q3D design with stacking and secondary loss models

- it does not include the database generation tools or the iterative solver

Many of these exist/are underway and will hopefully be included in subsequent versions.





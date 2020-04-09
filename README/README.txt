See also the accompanying README.ipynb, a Jupyter notebook with significantly more detailed explanation of the model.

The aim of this code is to model the aerodynamic performance of the helium turbine used in Reaction Engines' SABRE, which operates in an unconventional design space. While efforts have been made to generalise prediction ability, it should be noted that the methods applied here may not be compatible with significantly different specifictions, particularly the methods associated with gas properties, loss mechanisms and thermo-mechanical behaviour.

In general, the program uses basic turbomachinery principles to take design parameters and convert them to a complete turbine geometry, from which performance can be predicted with the use of correlations for loss mechanisms. It is possible to vary stage parameeters through the machine using splines and a number of control points, which can then be fed into optimisation routines to seek optimal performance. Machine is predicted by considering stress constraints, and resulting geometries can then be fed into further mechanical modelling, see ZK's work.

The program makes significant use of the numpy, scipy and math libraries, and also uses cvs and matplotlib for some optional functionalities. The provided folder should contain the following files:

- turbine.py - contains all of the functions for turbine geometry and performance prediction aside from those associated with blade profiles. Also contains a function for optimisation of a turbine design.

- prof_gen.py - contains functions for the generation of blade profiles written by CJC. Slight modifications have been made by WD, as well as the addition of a function for calculating blade section area and suction surface length and a function for tabulating these in .csv files.

- GUI.py - functions for plotting the blades and meridional views of turbines. Also contains a function that generates a GUI with matplotlib that can be used for rapid visualisation and testing of turbine designs. This GUI is quite slow and should be ported to a library more suited to GUI design.

- Profile_areas.csv and ss_grid.csv - tabulated data from the prof_gen functions, used when speed is necessary as interpolatiomn from these is faster than generating new sections every time.

- Data Mining.py, Plots.py - both conatin none of the prediction code, but have been used to investigate behaviours. Data Mining is where the functions from turbine.py were generally run and contains code to set inputs and print outputs and for extracting data for things like mesh generation with Autogrid. Also contains script for various optimisation routines and a sensitivity analysis. Plots as it suggests is where scripts for plotting relevant results have been kept. Both of these files are extremely messy and can be replaced with no loss in functionality.

If you have any more questions about what something does email wd269@cam.ac.uk. If any of the code looks to be formatted in a strange way it has likely been done so to conform with the python style guide set out in PEP 8.

Will
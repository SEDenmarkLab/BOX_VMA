Example ASO calculation in full_lib_with_commercial, run from command line in mollienv

1. molli grid --mlib ../in_silico_libraries/main_conformers/out1/conformers_no_linker.mlib

calculates bounding grid for conformer distributions

2. run alt_aso.py

calculates aso, filters out redundant conformers




The ASO/AEIF descriptors for modelling and the focused analogue descriptors for those workflows are also provided here. These were calculated with ccheminfolib. Modelling descriptors with zero variance or with correlation coefficient > 0.95 have been omitted.

In fa_lib_only, there are example scripts for calculating the kmeans (Kmeans.py) and generating the elbow plot (make_elbow.py).
PART 1: In Silico Library Generation and Preprocessing for Conformer Generation

The main library contains everything except for the focused analogues.

1. Run python generate_library.py

Parses chemdraws, makes combinatorial library, and minimizes the library with obabel UFF force field. Output has CBr2 linker. Saved raw molli output in in_silico_library folder, and UFF minimized structures in in_silico_library_uff_min. There are 96160 in silico ligands.

2. Run python move_stereochem_checker.py

Checks files for the correct chirality and move good files to in_silico_library_uff_min_checked1 folder. Some files fail RDKit's CIPlabelling for runtime error. These structures are huge and RDKit struggles with CIP assignments for the adamantyl substituents. The failed files (47 files) were excluded, and are found in in_silico_library_uff_min_checked1 folder. Only 3 files actually had problematic stereochemistry, which were not moved from in_silico_library_uff_min folder.

The total size of the checked library is 96110.
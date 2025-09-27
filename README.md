# Optimal-Kappa-and-Omega-Diagonal-Preconditioning
This repository contains the code and data associated with the the paper ``New Insights and Algorithms for Optimal Diagonal Preconditioning'' by Saeed Ghadimi, Woosuk L. Jung, Arnesh Sujanani, David Torregrosa-Belén, Henry Wolkowicz

The ``OurCodeVsStanford_Used_Paper'' folder contains the solver, table, and run files that appear in Sections 3.1 and 3.2 of the paper.
Please run 
1. table_Suite_Stan_three_comp.m to generate Table D.1
2. table_Rand_Stan_three_comp.m to generate Table D.2
3. table_Suite_Large_two_comp.m to generate Table D.3
4. table_Rand_Large_two_comp.m to generate Table D.4
5. table_StanvsUsPCGIn.m to generate Table D.5
6. table_PCGKappaSubgrad.m to generate Table D.6


The ''Arnesh-PCG'' folder contains the table files that eppear in Section 3.3 of the paper. Please run tablepaperKappaLinuxLarge.m and tablepaperKappaMedium.m
to generate the tables.

Also note that the data folder inside ''Arnesh-PCG'' folder contains .mat files of the SuiteSparse matrices that we consider in our experiments.


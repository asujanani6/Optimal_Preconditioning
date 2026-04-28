
As of April 22, 2026
This folder contains the table, tex, and run files for Tables 9-10 in the second draft (revision of the paper). These tables contain the experiments for PCG comparison between optimal kappa and optimal omega.

For example, to reproduce the table_9_PCG.tex, the user just needs to run table_9_PCG.m. Likewise, to reproduce table_10_PCG_Linux.tex, the user just needs to run table_10_PCG_Linux.m.






-----------------------

----------------------------------------------------
The ''Arnesh-PCG'' folder contains the table files that eppear in Section 3.3 of
the paper. Please run tablepaperKappaLinuxLarge.m and tablepaperKappaMedium.m to
generate the tables.
----------------------------------------------------
reorganizing/renaming files  table files for starting problem
tablepaper...m files for the tables in the paper
             large using linux server; medium for testing/changing
textable....tex files for the paper
Jun21/25
files reorganized; table files renamed e.g. tableLarge....
profile if statement used along with ispc to avoid linux errors
Jun17/25
table_testsOptKappa_Sparse_Large.m table_testsOptOmega_Sparse_Large.m
added for large/sparse problems to run on linux servers
also files polished and folder renamed/moved cleaned up
Jun13/25
Two matlab codes currently for two tables
opt kappa vs Jacobi
Jacobi with descent for kappa
table and run files using the pcg solver

Jun 16/25
for two large lables for Linux Servers, added two table files table_testsOptKappaLarge.m and table_testsOptOmegaLarge.m

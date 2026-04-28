# Optimal-Kappa-and-Omega-Diagonal-Preconditioning
This repository contains the code and data associated with the the paper ``Optimal Diagonal Preconditioning Beyond Worst-Case Conditioning: Theory and Practice of Omega Scaling'' by Saeed Ghadimi, Woosuk L. Jung, Arnesh Sujanani, David Torregrosa-Belén, Henry Wolkowicz

There are 3 important folders which contain the table and run files needed to reproduce the 10 tables that appear in the second draft of the paper.

1. Folder "Tables_1-5_Paper_One_Sided" contains the table, tex, and run files for Tables 1-5 in the second draft (revision of the paper). The subfolder "New_Tables_1_2_Draft2" contains the table and run files for Tables 1-2. The table and run files for Tables 3-5 are contained directly in the folder.  For example, to reproduce the table_3_Suite.tex, the user just needs to run table_3_Suite.m. Likewise, to reproduce table_4_Rand.tex, the user just needs to run table_4_Rand.m etc. The table and run files that are not used in the second draft of the paper but were used in the first draft of the paper are contained in the subfolders "Old_PCG_Tables_Draft_1" and "Old_Tables_1_2_Draft1".


2. Folder "Tables_6-8_Paper_Two_Sided" contains the table, tex, and run files for Tables 6-8 in the second draft (revision of the paper). These tables contain the experiments for the two-sided optimal omega vs optimal kappa comparison. For example, to reproduce the table_6_Suite.tex, the user just needs to run table_6_Suite.m. Likewise, to reproduce table_7_Rand.tex, the user just needs to run table_7_Rand.m etc.

3. Folder "Tables_9-10_Paper_PCG" contains the table, tex, and run files for Tables 9-10 in the second draft (revision of the paper). These tables contain the experiments for PCG comparison between optimal kappa and optimal omega. For example, to reproduce the table_9_PCG.tex, the user just needs to run table_9_PCG.m. Likewise, to reproduce table_10_PCG_Linux.tex, the user just needs to run table_10_PCG_Linux.m.


Note that each of these 3 folders also has its own Readme file containing the same information as above.

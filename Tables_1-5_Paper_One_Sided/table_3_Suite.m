%%%%%%%%%Table File for Minimizing Kappa on Large Instances: Comparison of
%%%%%%%%%Our Two Codes
%%%%%%%%

clc
clear 
close all

addpath(genpath('.'))

starttable = tic;
profilechoice = false;
if profilechoice && ispc
	profile clear
	profile on
end
fprintf('\nStarting Our Two Subgrad Methods\n')


%%%%%%%%%%%%%%%%%%%%%%
datavec= ["Pres_Poisson.mat", "bcsstk25.mat", "bcsstm25.mat", "gyro_m.mat", "gyro.mat", "bcsstk36.mat", "wathen100.mat", "wathen120.mat", "minsurfo.mat", "gridgena.mat", "G2_circuit.mat"];
%datavec= ["Pres_Poisson.mat", "bcsstk25.mat", "gyro_m.mat", "gyro.mat", "bcsstk36.mat", "wathen100.mat", "wathen120.mat", "minsurfo.mat", "gridgena.mat", "G2_circuit.mat"];


%%Filename for Table File
filename = 'table_3_Suite.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 


paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;


%%Call Run File
run_for_table_3(datavec,paramsUs,optionsUs,paramsUsSimp,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

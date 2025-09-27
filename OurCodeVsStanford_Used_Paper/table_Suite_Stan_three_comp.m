%%%%%%%%%Table File for Suite Sparse Stan vs Us vs Us Simplex
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
fprintf('\nStarting SuiteSparse Stan vs Our Two Subgrad Methods\n')


%%%%%%%%%%%%%%%%%%%%%%
datavec=["bcsstk08.mat", "bcsstk13.mat", "bcsstk21.mat", "bcsstk23.mat", "bcsstk24.mat", "bcsstk26.mat", "bcsstk28.mat", "bcsstk34.mat", "494_bus.mat", "662_bus.mat", "nasa1824.mat", "nasa2146.mat", "nasa2910.mat", "nos1.mat", "nos2.mat", "nos4.mat", "nos5.mat", "nos7.mat"];


%%Filename for Table File
filename = 'table_Suite_Stan_three_comp.tex';


paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 

%%Call Run File
run_Suite_Stan_three_comp(datavec,paramsUs,optionsUs,paramsStan,paramsUsSimp,filename)

if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

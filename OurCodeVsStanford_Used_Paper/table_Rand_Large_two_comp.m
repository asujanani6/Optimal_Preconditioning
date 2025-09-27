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
fprintf('\nStarting Our Subgrad Method\n')


%%%%%%%%%%%%%%%%%%%%%%
% dimvec=25000:15000:160000;
%dimvec=[50000,60000,70000,80000,100000,110000,120000,130000,140000,150000];
dimvec=50000:10000:150000;
%dimvec=20000:10000:40000;
seedvec=1:length(dimvec);

%%Filename for Table File
filename = 'table_Rand_Large_two_comp.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 


paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;


%%Call Run File
run_Rand_Large_two_comp(dimvec,seedvec,paramsUs,optionsUs,paramsUsSimp,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

%%%%%%%%%Table File for Opt Kappa vs Jacobi (Sparse Case)
%%% construct an A that is kappa diag-scaling optimal
%%% form J using A that is omega diag-scaling optimal
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
fprintf('\nStarting Random Stan vs Our Two Subgrad Method\n')


%%%%%%%%%%%%%%%%%%%%%%
dimvec=2000:500:15000;
%dimvec=1000:500:15000;
%dimvec=750:750:15000;
%dimvec=5000:500:6000;
seedvec=1:length(dimvec);
seedvec=seedvec*5;

%%Filename for Table File
filename = 'table_Rand_Stan_three_comp.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 

%%Call Run File
run_Rand_Stan_three_comp(dimvec,seedvec,paramsUs,optionsUs,paramsStan,paramsUsSimp,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

%%%%%%%%%Table File: Compare Random Instances Minimizing Kappa Among the 4
%%%%%%%%%Codes
%%%%%%%%
clc
clear 
close all

addpath(genpath('..'))
addpath /Users/arneshsujanani/Documents/mosek/11.1/toolbox/r2022b 
addpath('/Users/arneshsujanani/Documents/mosek/11.1/tools/platform/osxaarch64/bin')
addpath /Users/arneshsujanani/mosek


starttable = tic;
profilechoice = false;
if profilechoice && ispc
	profile clear
	profile on
end
fprintf('\nStarting Random Stan vs Our Two Subgrad Method\n')


%%%%%%%%%%%%%%%%%%%%%%
dimvec=2000:500:15000;
% dimvec=2000:500:5000;
seedvec=1:length(dimvec);
seedvec=seedvec*5;

%%Filename for Table File
filename = 'table_2_Rand.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 

paramsSub.ptype='S';
paramsSub.method="SDP";

%%%
%%Call Run File
run_for_table_2(dimvec,seedvec,paramsUs,optionsUs,paramsStan,paramsUsSimp,paramsSub,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

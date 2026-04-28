%%%%%%%%%Table File for PCG Comparsion: Us vs Stanford: Random
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
fprintf('\nStarting PCG Our Subgrad Method vs No Preconditioning\n')


%%%%%%%%%%%%%%%%%%%%%%
% dimvec=1000:300:15000;
dimvec=1000:300:10000;
%dimvec=1000:1000:5000;
seedvec=1:length(dimvec);

%%Filename for Table File
filename = 'tablePCGKappaSubgradA.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

%%Call Run File
run_PCGKappaSubgrad(dimvec,seedvec,paramsUs,optionsUs,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

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
fprintf('\nStarting Random PCG Stan vs Our Subgrad Method\n')


%%%%%%%%%%%%%%%%%%%%%%
%dimvec=1000:500:3000;
dimvec=1000:300:15000;
seedvec=1:length(dimvec);

%%Filename for Table File
filename = 'table_StanvsUsPCGIn.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 

%%Call Run File
run_StanvsUsPCGIn(dimvec,seedvec,paramsUs,optionsUs, paramsStan, filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

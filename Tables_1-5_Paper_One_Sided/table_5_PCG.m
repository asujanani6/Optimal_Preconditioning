%%%%%%%%%Table 5 PCG File
clear 
close all

addpath(genpath('.'))

%%%%%%%%%Check if CVX is available.
if ~exist('det_rootn')
	cvx_setup
end

starttable = tic;
profilechoice = false;
if profilechoice && ispc
	profile clear
	profile on
end
fprintf('\nStarting table_5_PCG\n')

%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
dimvec=1000:300:15000;
seedvec=1:length(dimvec);

%%Filename for Table File
filename = 'table_5_PCG.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 

%%Call Run File
run_for_table_5(dimvec,seedvec,paramsUs,optionsUs, paramsStan, filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('\nEnding table_5_PCG\n')
fprintf('total time table file %g\n',endtable)

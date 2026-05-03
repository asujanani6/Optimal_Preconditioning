%%%%%%%%%Table 2 Rand File
clear 
close all

addpath(genpath('..'))

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
fprintf('\nStarting table_2_Rand\n')


%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
dimvec=2000:500:15000;
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

%%Call Run File
run_for_table_2(dimvec,seedvec,paramsUs,optionsUs,paramsStan,paramsUsSimp,paramsSub,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('\nEnding table_2_Rand\n')
fprintf('total time table file %g\n',endtable)

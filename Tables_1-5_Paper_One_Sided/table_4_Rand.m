%%%%%%%%%Table 4 Rand File
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
fprintf('\nStarting table_4_Rand\n')

%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
dimvec=50000:10000:150000;
seedvec=1:length(dimvec);

%%Filename for Table File
filename = 'table_4_Rand.tex';

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 


paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;


%%Call Run File
run_for_table_4(dimvec,seedvec,paramsUs,optionsUs,paramsUsSimp,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('\nEnding table_4_Rand\n')
fprintf('total time table file %g\n',endtable)

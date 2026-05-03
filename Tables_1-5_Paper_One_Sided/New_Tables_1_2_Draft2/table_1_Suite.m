%%%%%%%%%Table 1 Suite File
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
fprintf('\nStarting table_1_Suite\n')

%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
datavec=["bcsstk08.mat", "bcsstk13.mat", "bcsstk21.mat", "bcsstk23.mat", "bcsstk24.mat", "bcsstk26.mat", "bcsstk28.mat", "bcsstk34.mat", "494_bus.mat", "662_bus.mat", "nasa1824.mat", "nasa2146.mat", "nasa2910.mat", "nos1.mat", "nos2.mat", "nos4.mat", "nos5.mat", "nos7.mat"];

%%Filename for Table File
filename = 'table_1_Suite.tex';

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
run_for_table_1(datavec,paramsUs,optionsUs,paramsStan,paramsUsSimp, paramsSub, filename)

if profilechoice && ispc
    profile report
end
endtable = toc(starttable);
fprintf('\nEnding table_1_Suite\n')
fprintf('total time table file %g\n',endtable)

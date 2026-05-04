%%%%%%%%%Table 3 Suite File
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

%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
datavec= ["Pres_Poisson.mat", "bcsstk25.mat", "bcsstm25.mat", "gyro_m.mat", "gyro.mat", "bcsstk36.mat", "wathen100.mat", "wathen120.mat", "minsurfo.mat", "gridgena.mat", "G2_circuit.mat"];

%%Filename for Table File
filename = 'table_3_Suite.tex';
startdatetime = datetime;
fprintf('\nStarting %s at %s\n',filename,startdatetime);

paramsUs.tolerance=1e-4;
paramsUs.maxitermain = 80;

optionsUs.plots = false;

paramsStan.ptype = 'R';
paramsStan.perturb = true; 


paramsUsSimp.hatdelta=1e-3;
paramsUsSimp.tolerance=1e-4;
paramsUsSimp.maxitermain=500;


%%Call Run File
run_for_table_3(datavec,paramsUs,optionsUs,paramsUsSimp,filename)


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)
enddatetime = datetime;
fprintf('\nEnding %s at %s\n',filename,enddatetime)
fprintf('\nStarting %s was at %s and ending at   %s\n', ...
	    filename,startdatetime,enddatetime)

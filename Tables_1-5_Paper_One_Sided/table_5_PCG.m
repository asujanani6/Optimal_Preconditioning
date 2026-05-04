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

%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
dimvec=1000:300:15000;
seedvec=1:length(dimvec);

%%Filename for Table File
filename = 'table_5_PCG.tex';
startdatetime = datetime;
fprintf('\nStarting %s at %s\n',filename,startdatetime);


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
fprintf('total time table file %g\n',endtable)
enddatetime = datetime;
fprintf('\nEnding %s at %s\n',filename,enddatetime)
fprintf('\nStarting %s was at %s and ending at   %s\n', ...
	    filename,startdatetime,enddatetime)

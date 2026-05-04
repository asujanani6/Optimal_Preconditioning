%%%%%%%%%Table 9 PCG File
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


%%%% SetTolerance for PCG
params.tol = 1e-7;
params.maxiter = 1e4;

%%%%%Load/Generate Problem Instances
dimvec=5000:5000:50000;
seedvec=1:length(dimvec);
densityvec = linspace(1e-1,1e-2,length(dimvec));
densityvec = max(densityvec,5./dimvec);
densityvec(end)=densityvec(end-1);

%%Specify how many initial points want to try/run for each problem instance
numofInitialPoints=5;



%%Filename for Table File
filename = 'table_9_PCG.tex';
startdatetime = datetime;
fprintf('\nStarting %s at %s\n',filename,startdatetime);

%%Call Run File
run_testsOptKappa_Sparse_PCG(dimvec,seedvec,densityvec,numofInitialPoints,params,filename);


if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)
enddatetime = datetime;
fprintf('\nEnding %s at %s\n',filename,enddatetime)
fprintf('\nStarting %s was at %s and ending at   %s\n', ...
	    filename,startdatetime,enddatetime)
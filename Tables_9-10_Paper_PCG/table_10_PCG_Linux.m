%%%%%%%%%Table 10 PCG Linux File
%%% construct an A that is kappa diag-scaling optimal
%%% form J using A that is omega diag-scaling optimal
%%%       to be run on fastlinux server
%cpu155.math.private 	Dell PowerEdge R650 	Two Intel Xeon Gold 6334 8-core 3.6 GHz (Ice Lake) 	256 GB
%%%%%%%%

clear 
close all

addpath(genpath('.'))

%%%%%%%%%Check if CVX is available.
if ~exist('det_rootn')
	cvx_setup
end

starttable = tic;
profilechoice = true;
if profilechoice && ispc
	profile clear
	profile on
end
fprintf('\nStarting table_10_PCG_Linux\n')


%%%% SetTolerance for PCG
params.tol = 1e-7;
params.maxiter = 1e4;

%%%%Load/Generate Problem Instances
%%Specify Seeds and densities (for sparse Opt Kappa A)
if ispc
	fprintf('Since a PC, SMALL size problems are being used\n')
	dimvec=500:20:900;
	seedvec=1:length(dimvec);
	exps = linspace(-1,1e-4,length(dimvec));
	densityvec = 10.^exps;
	densityvec = linspace(1e0,1e-3,length(dimvec));
	densityvec = max(densityvec,3./dimvec);
else
	fprintf('Since NOT a PC, linux exptected; LARGE size problems\n')
	dimvec=50000:2000:90000;
	seedvec=1:length(dimvec);
	densityvec = linspace(1e-1,5*1e-4,length(dimvec));
	densityvec = max(densityvec,3./dimvec);
end




%%Specify how many initial points want to try/run for each problem instance
numofInitialPoints=5;





%%Filename for Table File
filename = 'table_10_PCG_Linux.tex';

%%Call Run File
%[meanPCG_ratio_iter, meanPCG_ratio_time, densityA, condDiff, omegaDiff, dimA]= run_testsOptKappa_Sparse_PCG(dimvec,seedvec,densityvec,numofInitialPoints,params,filename);
run_testsOptKappa_Sparse_PCG(dimvec,seedvec,densityvec,numofInitialPoints,params,filename);
if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('\nEnding table_10_PCG_Linux\n')
fprintf('total time table file %g\n',endtable)

%%%%%%%%%Table File for Opt Kappa vs Jacobi (Sparse Case)
%%% construct an A that is kappa diag-scaling optimal
%%% form J using A that is omega diag-scaling optimal
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
fprintf('\nStarting OptKappaSparse tests\n')


%%%% SetTolerance for PCG
params.tol = 1e-7;
params.maxiter = 1e4;

%%Specify Seeds and densities (for sparse Opt Kappa A)
%densityvec=[1e-1, 5e-2, 1e-2 1e-3];  
%dimvec=2500:200:3100;  % for medium for testing% for medium testing
%seedvec=1:length(dimvec);
%tempd = 1:length(dimvec);
%densityvec = (1:length(dimvec)).*(10.^(-tempd));


% dimvec=500:20:700;
dimvec=5000:5000:50000;
seedvec=1:length(dimvec);
% densityvec = linspace(1e-1,1e-2,length(dimvec));
% densityvec = linspace(1e-1,1e-4,length(dimvec));
densityvec = linspace(1e-1,1e-2,length(dimvec));
% keyboard
densityvec = max(densityvec,5./dimvec);
densityvec(end)=densityvec(end-1);
% densityvec = max(densityvec,5./dimvec);

%seedvec=[1];
% densityvec=[1e-4, 1e-5, 1e-6];
% dimvec=[100000,150000,170000];  %%%%for large for testing large
% instances, uncomment for linux machines


%%Specify how many initial points want to try/run for each problem instance
numofInitialPoints=5;





%%Filename for Table File
filename = 'table_9_PCG.tex';

%%Call Run File
%[meanPCG_ratio_iter, meanPCG_ratio_time, densityA, condDiff, omegaDiff, dimA]= run_testsOptKappa_Sparse_PCG(dimvec,seedvec,densityvec,numofInitialPoints,params,filename);
run_testsOptKappa_Sparse_PCG(dimvec,seedvec,densityvec,numofInitialPoints,params,filename);
if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)

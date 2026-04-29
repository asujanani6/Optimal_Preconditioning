%%%%%%%%%Table File for Suite Sparse Stan vs Us vs Us Simplex
clc
clear
close all

addpath(genpath('..'))

% addpath /Users/arneshsujanani/Documents/mosek/11.1/toolbox/r2022b
% addpath('/Users/arneshsujanani/Documents/mosek/11.1/tools/platform/osxaarch64/bin')
% addpath /Users/arneshsujanani/mosek

%% Check that CVX is available
if exist('cvx_begin','file') ~= 2
    error(['CVX is not on the MATLAB path. ', ...
        'Please install CVX and run cvx_setup once before running this file.']);
end

%% Check that MOSEK is available through CVX
try
    cvx_solver mosek
catch ME
    error(['CVX is installed, but MOSEK is not available through CVX. ', ...
        'Run the following once in MATLAB: ', ...
        'cvx_setup; cvx_solver mosek; cvx_save_prefs. ', ...
        'Original error: %s'], ME.message);
end


starttable = tic;
profilechoice = false;
if profilechoice && ispc
    profile clear
    profile on
end
fprintf('\nStarting SuiteSparse Stan vs Our Two Subgrad Methods\n')


%%%%%%%%%%%%%%%%%%%%%%
datavec=["bcsstk08.mat", "bcsstk13.mat", "bcsstk21.mat", "bcsstk23.mat", "bcsstk24.mat", "bcsstk26.mat", "bcsstk28.mat", "bcsstk34.mat", "494_bus.mat", "662_bus.mat", "nasa1824.mat", "nasa2146.mat", "nasa2910.mat", "nos1.mat", "nos2.mat", "nos4.mat", "nos5.mat", "nos7.mat"];
% datavec=["bcsstk34.mat"];
% % datavec=["bcsstk08.mat"];
% datavec=["494_bus.mat"];

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
fprintf('total time table file %g\n',endtable)

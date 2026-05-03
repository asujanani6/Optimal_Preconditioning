%%%%%%%%%Table 8 Suite File
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
fprintf('\nStarting table_8_Suite\n')


%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
datavec=["airfoil_2d.mat", "c-24.mat", "c-36.mat", "c-39.mat", "cavity21.mat", "conf5_0-4x4-18.mat", "coupled.mat", "delaunay_n11.mat", "delaunay_n13.mat", "epb1.mat", "ex24.mat", "G34.mat", "G36.mat", "G59.mat", "G63.mat", "G64.mat", "garon1.mat", "Hamrle2.mat", "jan99jac040sc.mat", "Na5.mat", "t2d_q9.mat", "tols4000.mat", "utm5940.mat", "viscoplastic1.mat", "whitaker3.mat"];
seedvec=1:length(datavec);

%%%Our Solver Options
optionsOurSolver.maxiterscale = 2000;  % for left-right scaling steps while loop
optionsOurSolver.tolcd = 1e-2;        % stopping tol for scaling steps while loop
optionsOurSolver.verbose=false;


%%LSQR Options
optionsLSQR.maxitlsqr=100000;
optionsLSQR.tolsqr=1e-6;

%%%%Filename for Table File
filename='table_8_Suite.tex';

run_for_table_8(datavec,seedvec,optionsOurSolver,optionsLSQR,filename)

if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('\nEnding table_8_Suite\n')
fprintf('total time table file %g\n',endtable)

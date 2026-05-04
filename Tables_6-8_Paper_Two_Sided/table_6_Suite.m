%%%%%%%%%Table 6 Suite File
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

%%%Load/Generate Problem Instances
%%%%%%%%%%%%%%%%%%%%%%
datavec= ["08blocks.mat", "arc130.mat", "bwm200.mat", "ck104.mat", "ex1.mat", "football.mat", "gre_185.mat", "gre_216a.mat", "impcol_a.mat", "impcol_c.mat", "impcol_e.mat", "jazz.mat", "lshp_265.mat", "olm100.mat", "polbooks.mat", "rajat11.mat", "robot.mat", "rotor1.mat", "rw136.mat", "tub100.mat", "utm300.mat", "west0132.mat"];
seedvec=1:length(datavec);


%%%Our Solver Options
optionsOurSolver.maxiterscale = 2000;  % for left-right scaling steps while loop
optionsOurSolver.tolcd = 1e-2;        % stopping tol for scaling steps while loop
optionsOurSolver.verbose=false;


%%%%Stanford Parameter Options
paramsStan.ptype = 'T';

%%LSQR Options
optionsLSQR.maxitlsqr=5000;
optionsLSQR.tolsqr=1e-8;

%%Filename for Table File
filename='table_6_Suite.tex';
startdatetime = datetime;
fprintf('\nStarting %s at %s\n',filename,startdatetime);

run_for_table_6(datavec,seedvec,optionsOurSolver,paramsStan,optionsLSQR,filename)

if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)
enddatetime = datetime;
fprintf('\nEnding %s at %s\n',filename,enddatetime)
fprintf('\nStarting %s was at %s and ending at   %s\n', ...
	    filename,startdatetime,enddatetime)

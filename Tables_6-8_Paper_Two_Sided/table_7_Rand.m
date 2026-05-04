%%%%%%%%%Table 7 Rand File
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
%%%Dim, Seed, Density, Reciprocal Condition Number
dimvec=linspace(50, 300, 26);
seedvec=1:length(dimvec);

for i=1:length(dimvec)
densityvec(i) = 1/(log(dimvec(i)));
end
rcvec=linspace(1e-3, 1e-5, 26);

%%%Our Solver Options
optionsOurSolver.maxiterscale = 2000;  % for left-right scaling steps while loop
optionsOurSolver.tolcd = 1e-2;        % stopping tol for scaling steps while loop
optionsOurSolver.verbose=false;


%%%%Stanford Parameter Options
paramsStan.ptype = 'T';

%%LSQR Options
optionsLSQR.maxitlsqr=5000;
optionsLSQR.tolsqr=1e-8;

%%%%
filename='table_7_Rand.tex';
startdatetime = datetime;
fprintf('\nStarting %s at %s\n',filename,startdatetime);

run_for_table_7(dimvec,seedvec,densityvec,rcvec,optionsOurSolver,paramsStan,optionsLSQR,filename)

if profilechoice && ispc
	profile report
end
endtable = toc(starttable);
fprintf('total time table file %g\n',endtable)
enddatetime = datetime;
fprintf('\nEnding %s at %s\n',filename,enddatetime)
fprintf('\nStarting %s was at %s and ending at   %s\n', ...
	    filename,startdatetime,enddatetime)

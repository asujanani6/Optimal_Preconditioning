%%%%Table File to compare Stanford vs US for Two-Sided LSQR Comp


%%%Dim, Seed, Density, Reciprocal Condition Number

datavec= ["08blocks.mat", "arc130.mat", "bwm200.mat", "ck104.mat", "ex1.mat", "football.mat", "gre_185.mat", "gre_216a.mat", "impcol_a.mat", "impcol_c.mat", "impcol_e.mat", "jazz.mat", "lshp_265.mat", "olm100.mat", "polbooks.mat", "rajat11.mat", "robot.mat", "rotor1.mat", "rw136.mat", "tub100.mat", "utm300.mat", "west0132.mat"];
% datavec= ["ex1.mat"];
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

%%%%
filename='table_6_Suite.tex'
% filename='singlematrix.tex';


run_for_table_6(datavec,seedvec,optionsOurSolver,paramsStan,optionsLSQR,filename)
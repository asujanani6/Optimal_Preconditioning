%%%%Table File to compare Stanford vs US for Two-Sided LSQR Comp


%%%Dim, Seed, Density, Reciprocal Condition Number
% dimvec=linspace(100, 300, 21);
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
filename='table_7_Rand.tex'

run_for_table_7(dimvec,seedvec,densityvec,rcvec,optionsOurSolver,paramsStan,optionsLSQR,filename)

function run_for_table_8(datavec,seedvec,optionsOurSolver,optionsLSQR, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Two-Sided Comparison between our code and stanford group's code

%%%Number of Problems is just number of dimensions inputted.
numofProblems=length(datavec);

%%%Save Densities of A matrix
densityA=zeros(numofProblems,1);

%%%Save dimensions of A matrix
dimA=zeros(numofProblems,1);

%%%Save kappa and omega condition numbers of A
omegaA=zeros(numofProblems,1);



%%%Save CPU time: Us
US_prec_time=zeros(numofProblems,1);


%%%Save LSQR time: Us vs No Prec
US_LSQR_time=zeros(numofProblems,1);
None_LSQR_time=zeros(numofProblems,1);

%%%Total CPU time: Us
US_tot_time=zeros(numofProblems,1);


%%%%% vectors to save output in terms of ratios of condition numbers
omega_USDividedA=zeros(numofProblems,1);


%%%Save LSQR iterations: Us vs A
US_LSQR_iter=zeros(numofProblems,1);
A_LSQR_iter=zeros(numofProblems,1);


for i=1:numofProblems
    startgen=tic;
    %%%%Generate/Get Problem
    load(datavec(i));
    A=Problem.A;
    n=size(A,2);


    rng(seedvec(i));
    b=randn(n,1);
    omegaa = @(M)(  (trace(M'*M)/n)/(det_rootn(M'*M))  );

    %%%%Save Dimensions, Densities,Omega
    dimA(i)=n;
    densityA(i)=nnz(A)/numel(A);

    endgen=toc(startgen);

    fprintf('\nNEW prob: time to gen = %g; size n =%i; density %g \n', ...
        endgen,dimA(i),densityA(i))


    omegaA(i)=omegaa(A);



    %%%%Run Our Code
    [csqs,dsqs,output] = Solver_Alt_TwoSided(A,optionsOurSolver);
    US_prec_time(i)=output.time;

    M1 = spdiags(csqs,0,n,n);
    M2 = spdiags(dsqs,0,n,n);
    M1AM2=M1\A/M2;

    %%%Compute Ratio of Omega Values
    omega_USDividedA(i)=omegaa(M1AM2)/omegaa(A);

    %%%%Solve LSQR: Compare Our Preconditioner vs Theirs
    tic
    [x2,flag2,relres2,iter2] = lsqr(( M1\A/M2) , M1\b, optionsLSQR.tolsqr, optionsLSQR.maxitlsqr, ...
        speye(n), speye(n));
    US_LSQR_time(i)=toc;

    tic
    [x3,flag3,relres3,iter3] = lsqr((A) , b, optionsLSQR.tolsqr, optionsLSQR.maxitlsqr, ...
        speye(n), speye(n));
    None_LSQR_time(i)=toc;

    US_LSQR_iter(i)=iter2;
    A_LSQR_iter(i)=iter3;

    US_tot_time(i)=US_LSQR_time(i)+US_prec_time(i);

end

%%%%Create LaTeX Table

fmt1 = '%5.0f & %6.1e & %6.1e & ';


fmt3 = '%6.1e & %6.1e & %6.0f & %6.0f & %5.2e & %5.2e';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-3}\\cline{4-9}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|ccc|cc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{3}{|c||}{Dim \& Density, \& $\omega(A)$} & ', ...
    '\multicolumn{2}{|c|}{Ratio of Omega \& Prec CPU Time} & ', ...
    '\multicolumn{2}{|c|}{LSQR Iter} & ', ...
    '\multicolumn{2}{|c|}{Total CPU Time}', ...
    ]);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density', '$\omega(A)$');
fprintf(fid, hfmt3,'$\omega(S)/\omega(A)$','\Cref{alg:Two-Sided} CPU', ...
    'No Prec', '\Cref{alg:Two-Sided}', 'No Prec', '\Cref{alg:Two-Sided}');
fprintf(fid, fmt4);

for k=1:numofProblems
    fprintf(fid, fmt1, dimA(k), densityA(k), omegaA(k));
    fprintf(fid, fmt3, omega_USDividedA(k), US_prec_time(k), A_LSQR_iter(k), US_LSQR_iter(k), None_LSQR_time(k), US_tot_time(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

type(filename)


end



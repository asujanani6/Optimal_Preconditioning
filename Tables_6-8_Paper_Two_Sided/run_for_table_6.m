function run_for_table_6(datavec,seedvec,optionsOurSolver,paramsStan,optionsLSQR, filename)

addpath(genpath('..'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Two-Sided Comparison between our code and stanford group's code

%%%Number of Problems is just number of dimensions inputted.
numofProblems=length(datavec);

%%%Save Densities of A matrix
densityA=zeros(numofProblems,1);

%%%Save dimensions of A matrix
dimA=zeros(numofProblems,1);

%%%Save kappa and omega condition numbers of A
kappaA=zeros(numofProblems,1);
omegaA=zeros(numofProblems,1);



%%%Save CPU time: Us vs Stanford
US_prec_time=zeros(numofProblems,1);
Stan_prec_time=zeros(numofProblems,1);

%%Save LSQR time: Us vs Stanford
US_LSQR_time=zeros(numofProblems,1);
Stan_LSQR_time=zeros(numofProblems,1);

%%%Save Total CPU time: Us vs Stanford
US_total_time=zeros(numofProblems,1);
Stan_total_time=zeros(numofProblems,1);

%%%%% vectors to save output in terms of ratios of condition numbers
omega_USDividedStan=zeros(numofProblems,1);
kappa_USDividedStan=zeros(numofProblems,1);



%%%Save LSQR iterations: Us vs Stanford
US_LSQR_iter=zeros(numofProblems,1);
Stan_LSQR_iter=zeros(numofProblems,1);


for i=1:numofProblems

%%%%Generate/Get Problem
 load(datavec(i));
 A=Problem.A;
 n=size(A,2);


rng(seedvec(i));
b=randn(n,1);
omegaa = @(M)(  (trace(M'*M)/n)/(det_rootn(M'*M))  );

%%%%Save Dimensions, Densities, Kappa, Omega
dimA(i)=n;
densityA(i)=nnz(A)/numel(A);

svdA=svd(full(A));
kappaA(i)=max(svdA)/min(svdA);
omegaA(i)=omegaa(A);



%%%%Run Our Code
[csqs,dsqs,output] = Solver_Alt_TwoSided(A,optionsOurSolver);
US_prec_time(i)=output.time;

M1 = spdiags(csqs,0,n,n);
M2 = spdiags(dsqs,0,n,n);
M1AM2=M1\A/M2;

%%%%Run Stanford's Code
Eprob = getoptprob(A, paramsStan);

% Solve problem using Stanford's code
tic
Eopt = optprecond(Eprob);
Stan_prec_time(i)=toc;


D=spdiags(Eopt.D,0,n,n);
E=spdiags(Eopt.E,0,n,n);
DAE=D*A*E;


%%%Compute Ratio of Omega Values
omega_USDividedStan(i)=omegaa(M1AM2)/omegaa(DAE);

%%%Compute Kappa Values and Ratio of Values
svdM1AM2=svd(full((M1AM2)));
svdDAE=svd(full((DAE)));
kappa_US_M1AM2=max(svdM1AM2)/min(svdM1AM2);
kappa_Stan_DAE=max(svdDAE)/min(svdDAE);

kappa_USDividedStan(i)=kappa_US_M1AM2/kappa_Stan_DAE;

%%%%Solve LSQR: Compare Our Preconditioner vs Theirs
tic
[x2,flag2,relres2,iter2] = lsqr(( M1\A/M2) , M1\b, optionsLSQR.tolsqr, optionsLSQR.maxitlsqr, ...
	          speye(n), speye(n));
US_LSQR_time(i)=toc;

tic
[x3,flag3,relres3,iter3] = lsqr(( D*A*E) , D*b, optionsLSQR.tolsqr, optionsLSQR.maxitlsqr, ...
	          speye(n), speye(n));
Stan_LSQR_time(i)=toc;


US_LSQR_iter(i)=iter2
Stan_LSQR_iter(i)=iter3

US_total_time(i)=US_prec_time(i)+US_LSQR_time(i);
Stan_total_time(i)=Stan_prec_time(i)+Stan_LSQR_time(i);

end

%%%%Create LaTeX Table

fmt1 = '%5.0f & %6.1e & %6.1e & ';


fmt3 = '%6.1e & %6.1e & %5.3f & %5.3f & %4.0f & %4.0f & %5.3f & %5.2f';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-3}\\cline{4-11}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|ccc|cc|cc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{3}{|c||}{Dim \& Density, \& $\kappa(A)$} & ', ...
    '\multicolumn{2}{|c|}{Ratio of Condition Numbers} & ', ...
    '\multicolumn{2}{|c|}{CPU Time for Prec} & ', ...
    '\multicolumn{2}{|c|}{LSQR Iterations} & ', ...
    '\multicolumn{2}{|c|}{Total CPU Time}'
    ]);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density', '$\kappa(A)$');
fprintf(fid, hfmt3,'$\kappa(S)/\kappa(T)$','$\omega(S)/\omega(T)$','\Cref{alg:Two-Sided}', ...
     '\cite{doi:10.1287/opre.2022.0592}', '\Cref{alg:Two-Sided}', '\cite{doi:10.1287/opre.2022.0592}', '\Cref{alg:Two-Sided}', '\cite{doi:10.1287/opre.2022.0592}');
fprintf(fid, fmt4);

for k=1:numofProblems
    fprintf(fid, fmt1, dimA(k), densityA(k), kappaA(k));
    fprintf(fid, fmt3, kappa_USDividedStan(k), omega_USDividedStan(k), US_prec_time(k), Stan_prec_time(k), US_LSQR_iter(k), Stan_LSQR_iter(k), US_total_time(k), Stan_total_time(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);


end




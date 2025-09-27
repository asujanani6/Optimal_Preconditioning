function run_Rand_Large_two_comp(dimvec,seedvec,paramsUs,optionsUs,paramsUsSimp,filename)

%%%%%%Minimize kappa comparsion between our code and stanford group's code

%%%Number of Problems is just number of dimensions inputted.
numofProblems=length(dimvec);


%%%%% vectors to save output. Reduction in Kappa: Us
kappa_Reduc_US=zeros(numofProblems,1);
kappa_Reduc_UsSimp=zeros(numofProblems,1);


%%%%% vectors to save output.  Kappa: Us
kappa_US=zeros(numofProblems,1);
kappa_UsSimp=zeros(numofProblems,1);


%%%Save CPU time: Us
US_kappa_time=zeros(numofProblems,1);
UsSimp_kappa_time=zeros(numofProblems,1);

%%%Save Densities of A matrix
densityA=zeros(numofProblems,1);

%%%Save dimensions of A matrix
dimA=zeros(numofProblems,1);


%%%Save kappa/condition number of A
kappaA=zeros(numofProblems,1);


for i=1:numofProblems

    %%%Loading Problem
    startgen = tic;

    rng(seedvec(i));

    n=dimvec(i);
    rc=5e-7/n;
    %rc=1e-8/n;
    % rc=1/n;
    density=(2/n)+(0.01/(n*log(n)));
    % density=0.01;
    A=sprandsym(n,density,rc,1);


    %%%%Save densities and dimensions of A
    dimA(i)=n;
    densityA(i)=nnz(A)/numel(A);

    endgen = toc(startgen);

    fprintf('NEW prob: time to gen = %g; size n =%i; density %g \n', ...
        endgen,n,densityA(i))

    R = chol(A);

    paramsUs.R = R;
    Afunctionalmax = @(x)( A*x );   % for max eig
    Afunctionalmin = @(x)( (R\(R'\x)) );   % for min eig
    [x1, lam1]= eigs(Afunctionalmax ,n,1,'LM', 'Tolerance', 1e-10);
    [xn, lamn] = eigs(Afunctionalmin, n, 1, 'SM', 'Tolerance', 1e-10);
    paramsUs.x1 = x1;
    paramsUs.xn = xn;

    kappaA(i) = lam1/lamn;
    steplgthmin = max(1e-14,min(1e-10/n,1/(n*kappaA(i))));

    paramsUs.steplgthmin = steplgthmin;

    tic
    output = OurSubgradSolver(A,paramsUs,optionsUs);
    US_kappa_time(i) = toc;

    w=output.w;

    Afunctionalmax = @(x)( A*(w.*x) );   % for max eig
    Afunctionalmin = @(x)(R\(R'\x))./w;  % for min eig


    [x1new, lam1new]= eigs(Afunctionalmax ,n,1, 'LM', 'StartVector', x1, 'Tolerance', 1e-10);

    [xnnew, lamnnew] = eigs(Afunctionalmin, n, 1, 'SM', 'StartVector', xn, 'Tolerance', 1e-10);

    kappa_US(i)=lam1new/lamnnew;

    kappa_Reduc_US(i)=((kappaA(i)-kappa_US(i))/kappaA(i))*100;


    paramsUsSimp.R=R;
    paramsUsSimp.x1 = x1;
    paramsUsSimp.xn = xn;

    tic
    outputSimplex = OurSubgradSolverSimplex(A,paramsUsSimp);
    UsSimp_kappa_time(i) = toc;



    wSimp=outputSimplex.w;

    AfunctionalmaxSimp = @(x)( A*(wSimp.*x) );   % for max eig
    AfunctionalminSimp = @(x)(R\(R'\x))./wSimp;   % for min eig


    [x1Simp, lam1Simp]= eigs(AfunctionalmaxSimp,n,1, 'LM', 'StartVector', x1, 'Tolerance', 1e-10);

    [xnSimp, lamnSimp] = eigs(AfunctionalminSimp, n, 1, 'SM', 'StartVector', xn, 'Tolerance', 1e-10);

    kappa_UsSimp(i)=lam1Simp/lamnSimp;

    kappa_Reduc_UsSimp(i)=((kappaA(i)-kappa_UsSimp(i))/kappaA(i))*100;


end




%%%%Create LaTeX Table

fmt1 = '%5.0f & %6.1e & %6.1e & ';


fmt3 = '%6.1e & %6.1e & %5.3f & %5.3f';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-3}\\cline{3-7}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|ccc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{3}{|c||}{Dim \& Density, \& $\kappa(A)$} & ', ...
    '\multicolumn{2}{|c|}{\% Reduction in $\kappa$} & ', ...
    '\multicolumn{2}{|c|}{CPU Time}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density', '$\kappa(A)$');
fprintf(fid, hfmt3,'Alg 2.2','Alg. 3.1','Alg. 2.2', ...
     'Alg 3.1');
fprintf(fid, fmt4);

for k=1:numofProblems
    fprintf(fid, fmt1, dimA(k), densityA(k), kappaA(k));
    fprintf(fid, fmt3, kappa_Reduc_UsSimp(k), kappa_Reduc_US(k), UsSimp_kappa_time(k), US_kappa_time(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);


end

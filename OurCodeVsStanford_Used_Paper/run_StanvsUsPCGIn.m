function run_StanvsUsPCGIn(dimvec,seedvec,paramsUs,optionsUs, paramsStan, filename)

%%%%%%Minimize kappa comparsion between our code and stanford group's code

%%%Number of Problems is just number of dimensions inputted.
numofProblems=length(dimvec);


%%%%% vectors to save output. Reduction in Kappa: Us vs Stanford
kappa_Reduc_US=zeros(numofProblems,1);
kappa_Reduc_Stan=zeros(numofProblems,1);


%%%%% vectors to save output.  Kappa: Us vs Stanford
kappa_US=zeros(numofProblems,1);
kappa_Stan=zeros(numofProblems,1);


%%%Save CPU time: Us vs Stanford
US_kappa_time=zeros(numofProblems,1);
Stan_kappa_time=zeros(numofProblems,1);


%%%Save PCG iterations: Us vs Stanford
US_PCG_iter=zeros(numofProblems,1);
Stan_PCG_iter=zeros(numofProblems,1);


%%%Save PCG CPU time: Us vs Stanford
US_PCG_CPU=zeros(numofProblems,1);
Stan_PCG_CPU=zeros(numofProblems,1);


%%Ratios

Ratio_kappa=zeros(numofProblems,1);
Ratio_CPU=zeros(numofProblems,1);



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
    % rc=1e-3/n;
    rc=1e-6/n;
    density=(2/n)+(10/(n*log(n)));
    %density=(2/n)+(1e-1/(n*log(n)));
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


    Eprob = getoptprob(R, paramsStan);

    % Solve problem using Stanford's code
    tic
    Eopt = optprecond(Eprob);
    Stan_kappa_time(i)=toc;


    EAE=(spdiags(Eopt.E, 0, n, n))*A*(spdiags(Eopt.E, 0, n, n));

    [x1EAE, lam1EAE]= eigs(EAE ,1, 'LM', 'Tolerance', 1e-10);
    [xnEAE, lamnEAE]= eigs(EAE ,1, 'SM', 'Tolerance', 1e-10);

    kappa_Stan(i)=lam1EAE/lamnEAE;



    kappa_Reduc_Stan(i)=((kappaA(i)-kappa_Stan(i))/kappaA(i))*100;

    Ratio_kappa(i)=kappa_Stan(i)/kappa_US(i);

    Ratio_CPU(i)=Stan_kappa_time(i)/US_kappa_time(i);

    b=randn(n,1);
    x0Start=randn(n,1);

    startUSpcg = tic;
    [xUS,flagUS,relresUS,iterUS]=pcg(A,b,1e-4,5e6,spdiags(1./w, 0, n, n),[],x0Start);
    US_PCG_CPU(i)=toc(startUSpcg);

    % flagUS
    % relresUS

    US_PCG_iter(i)=iterUS;


    startSTpcg = tic;
    [xST,flagST,relresST,iterST]=pcg(A,b,1e-4,5e6,spdiags(1./((Eopt.E).^2), 0, n, n),[],x0Start);
    Stan_PCG_CPU(i)=toc(startSTpcg);

    % flagST
    % relresST

    Stan_PCG_iter(i)=iterST;





    % tic
    % xPCGnew = pcg(A,b,1e-5,1e4,WPCG);
    % toc


end




%%%%Create LaTeX Table

fmt1 = '%5.0f & %6.1e & %6.1e & ';


fmt3 = '%6.1e & %6.1e & %5.0f & %5.0f & %4.2f & %4.2f';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-3}\\cline{3-9}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|ccc|cc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{3}{|c||}{Dim, Density, \& $\kappa(A)$} & ', ...
    '\multicolumn{2}{|c|}{\% Reduction in $\kappa$} & ', ...
    '\multicolumn{2}{|c|}{PCG Iterations} &', ...
    '\multicolumn{2}{|c|}{PCG CPU Time}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density', '$\kappa(A)$');
fprintf(fid, hfmt3,'Qu et al.', 'Alg. 3.1', ...
    'Qu et al.', 'Alg 3.1', ...
    'Qu et al.', 'Alg. 3.1');
fprintf(fid, fmt4);

for k=1:numofProblems
    fprintf(fid, fmt1, dimA(k), densityA(k), kappaA(k));
    fprintf(fid, fmt3, kappa_Reduc_Stan(k), kappa_Reduc_US(k), ...
        Stan_PCG_iter(k), US_PCG_iter(k), ...
        Stan_PCG_CPU(k), US_PCG_CPU(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);


end

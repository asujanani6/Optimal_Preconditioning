function run_for_table_1(datavec,paramsUs,optionsUs,paramsStan,paramsUsSimp, paramsSub, filename)

%%%%%%Minimize kappa comparsion between our code and stanford group's code

%%%Number of Problems is just number of dimensions inputted.
numofProblems=length(datavec);


%%%%% vectors to save output. Reduction in Kappa: Us vs Stanford
kappa_Reduc_US=zeros(numofProblems,1);
kappa_Reduc_UsSimp=zeros(numofProblems,1);
kappa_Reduc_Stan=zeros(numofProblems,1);
kappa_Reduc_Sub=zeros(numofProblems,1);


%%%%% vectors to save output.  Kappa: Us vs Stanford
kappa_US=zeros(numofProblems,1);
kappa_UsSimp=zeros(numofProblems,1);
kappa_Stan=zeros(numofProblems,1);
kappa_Sub=zeros(numofProblems,1);

%%%Save CPU time: Us vs Stanford
US_kappa_time=zeros(numofProblems,1);
UsSimp_kappa_time=zeros(numofProblems,1);
Stan_kappa_time=zeros(numofProblems,1);
Sub_kappa_time=zeros(numofProblems,1);

%%Ratios
Ratio_CPUSimp=zeros(numofProblems,1);
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

    load(datavec(i));
    A=Problem.A;
    n=size(A,2);

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
    [x1, lam1, flagA1]= eigs(Afunctionalmax ,n,1,'LM', 'Tolerance', 1e-10);
    [xn, lamn, flagAn] = eigs(Afunctionalmin, n, 1, 'SM','Tolerance', 1e-10);

    if flagA1==1
        [x1, lam1, flagA1]= eigs(Afunctionalmax ,n,1,'LM', 'Tolerance', 1e-6);
    end

    if flagAn==1
        [xn, lamn, flagAn] = eigs(Afunctionalmin, n, 1, 'SM','Tolerance', 1e-6);
    end

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


    [x1new, lam1new, flagnew1]= eigs(Afunctionalmax ,n,1, 'LM', 'StartVector', x1, 'Tolerance', 1e-10);

    [xnnew, lamnnew, flagnewn] = eigs(Afunctionalmin, n, 1, 'SM', 'StartVector', xn, 'Tolerance', 1e-10);

    if flagnew1==1
        [x1new, lam1new, flagnew1]= eigs(Afunctionalmax ,n,1, 'LM', 'StartVector', x1, 'Tolerance', 1e-6);
    end

    if flagnewn==1
        [xnnew, lamnnew, flagnewn] = eigs(Afunctionalmin, n, 1, 'SM', 'StartVector', xn, 'Tolerance', 1e-6);
    end

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


    [x1Simp, lam1Simp, flagsimp1]= eigs(AfunctionalmaxSimp,n,1, 'LM', 'StartVector', x1, 'Tolerance', 1e-10);

    [xnSimp, lamnSimp, flagsimpn] = eigs(AfunctionalminSimp, n, 1, 'SM', 'StartVector', xn, 'Tolerance', 1e-10);

    if flagsimp1==1
        [x1Simp, lam1Simp, flagsimp1]= eigs(AfunctionalmaxSimp,n,1, 'LM', 'StartVector', x1, 'Tolerance', 1e-6);
    end

    if flagsimpn==1
        [xnSimp, lamnSimp, flagsimpn] = eigs(AfunctionalminSimp, n, 1, 'SM', 'StartVector', xn, 'Tolerance', 1e-6);
    end

    kappa_UsSimp(i)=lam1Simp/lamnSimp;

    kappa_Reduc_UsSimp(i)=((kappaA(i)-kappa_UsSimp(i))/kappaA(i))*100;





    Eprob = getoptprob(R, paramsStan);

    % Solve problem using Stanford's code
    tic
    Eopt = optprecond(Eprob);
    Stan_kappa_time(i)=toc;



    sumNanVal=sum(isnan(Eopt.E));
    if sumNanVal>=1 || isreal(Eopt.pX)==0
        kappa_Stan(i)=kappaA(i);
    else
        EAE=(spdiags(Eopt.E, 0, n, n))*A*(spdiags(Eopt.E, 0, n, n));

        [x1EAE, lam1EAE, flagEAE1]= eigs(EAE ,1, 'LM', 'Tolerance', 1e-10);
        [xnEAE, lamnEAE, flagEAEn]= eigs(EAE ,1, 'SM', 'Tolerance', 1e-10);

        if flagEAE1==1
            [x1EAE, lam1EAE, flagEAE1]= eigs(EAE ,1, 'LM', 'Tolerance', 1e-6);
        end

        if flagEAEn==1
            [xnEAE, lamnEAE, flagEAEn]= eigs(EAE ,1, 'SM', 'Tolerance', 1e-6);
        end

        kappa_Stan(i)=lam1EAE/lamnEAE;

    end

    kappa_Reduc_Stan(i)=((kappaA(i)-kappa_Stan(i))/kappaA(i))*100;

    %%%%%%%%New for Revision: Subspace Code Testing
    paramsSub.subspace=[ones(n, 1), full(diag(A))];
    DSprob = getoptprob(A, paramsSub);

    tic
    DSopt = optprecondSUB(DSprob);
    Sub_kappa_time(i)=toc;


    sumNanValSub=sum(isnan(DSopt.E));
    if sumNanValSub>=1 || isreal(DSopt.pX)==0
        kappa_Sub(i)=kappaA(i);
    else
        DsADs=sqrt((spdiags(DSopt.E, 0, n, n)))*A*sqrt(spdiags(DSopt.E, 0, n, n));
        [x1DsADs, lam1DsADs, flagDs1]= eigs(DsADs ,1, 'LM', 'Tolerance', 1e-10);
        [xnDsADs, lamnDsADs, flagDsn]= eigs(DsADs ,1, 'SM', 'Tolerance', 1e-10);

        if flagDs1==1
            [x1DsADs, lam1DsADs, flagDs1]= eigs(DsADs ,1, 'LM', 'Tolerance', 1e-6);
        end

        if flagDsn==1
            [xnDsADs, lamnDsADs, flagDsn]= eigs(DsADs ,1, 'SM', 'Tolerance', 1e-6);
        end
        kappa_Sub(i)=lam1DsADs/lamnDsADs;

    end

    kappa_Reduc_Sub(i)=((kappaA(i)-kappa_Sub(i))/kappaA(i))*100;

    Ratio_CPUSimp(i)=Stan_kappa_time(i)/UsSimp_kappa_time(i);

    Ratio_CPU(i)=Stan_kappa_time(i)/US_kappa_time(i);

    % keyboard



end




%%%%Create LaTeX Table

fmt1 = '%5.0f & %6.1e & %6.1e & ';


fmt3 = '%6.1e & %6.1e & %6.1e & %6.1e & %5.3f & %5.3f & %5.3f & %5.3f & %4.1f & %4.1f';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-3}\\cline{4-13}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|ccc|cccc|cccc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{3}{|c||}{Dim \& Density, \& $\kappa(A)$} & ', ...
    '\multicolumn{4}{|c|}{\% Reduction in $\kappa$} & ', ...
    '\multicolumn{4}{|c|}{CPU Time (seconds)} &', ...
    '\multicolumn{2}{|c|}{CPU Ratios}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density', '$\kappa(A)$');
fprintf(fid, hfmt3, 'Qu et al.', 'Gao et al.', 'Alg 2.2', 'Alg. A.1', ...
    'Qu et al.', 'Gao et al.', 'Alg 2.2', 'Alg A.1', ...
    'Qu/Alg 2.1', 'Qu/Alg A.1');
fprintf(fid, fmt4);

for k=1:numofProblems
    fprintf(fid, fmt1, dimA(k), densityA(k), kappaA(k));
    fprintf(fid, fmt3, kappa_Reduc_Stan(k), kappa_Reduc_Sub(k), kappa_Reduc_UsSimp(k), kappa_Reduc_US(k), ...
        Stan_kappa_time(k), Sub_kappa_time(k), UsSimp_kappa_time(k), US_kappa_time(k), ...
        Ratio_CPUSimp(k),Ratio_CPU(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);


end

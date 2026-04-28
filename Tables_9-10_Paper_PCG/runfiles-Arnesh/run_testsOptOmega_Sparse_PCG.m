function run_testsOptOmega_Sparse_PCG(dimvec, seedvec, densityvec, rcvec, numofInitialPoints, paramsPCG, paramsGD, filename)
%%% Start with an omega-optimal scaled A and compare results
%%%   with kappa-scaled A obtained using steepest descent to reduce kappa


%%    det_rootn(B) is used below so cvx is needed
if ~exist('det_rootn')
	cvx_setup
end

%%%Number of Problems
numofProblems=length(dimvec);

%%Number of Trials (with Many Different Initial Points).

numofTrials=numofInitialPoints*numofProblems;

%%%%% vectors to save output. The difference
%%%%% between classical condition number of Jacobi and Jacobi descent is saved
condDiff=zeros(numofProblems,1);
condRatios=zeros(numofProblems,1);

%For each problem, the difference
%%%%% between omega condition number Jacobi and Jacobi descent is saved
omegaDiff=zeros(numofProblems,1);
omegaRatios=zeros(numofProblems,1);

%%%Save PCG iterations Jacobi and Jacobi+Kappa Descent as Well as Ratios
PCG_J_iter=zeros(numofTrials,1);
PCG_JK_iter=zeros(numofTrials,1);
PCG_ratio_iter=zeros(numofTrials,1);


%%%Save PCG time Jacobi and Jacobi+Kappa Descent as Well as Ratios
PCG_J_time=zeros(numofTrials,1);
PCG_JK_time=zeros(numofTrials,1);
PCG_ratio_time=zeros(numofTrials,1);


%%Save Mean Iterations over many different Initial Points for Each Problem
%%Instance
meanPCG_J_iter=zeros(numofProblems,1);
meanPCG_JK_iter=zeros(numofProblems,1);
meanPCG_ratio_iter=zeros(numofProblems,1);

%%Save Mean Time over many different Initial Points for Each Problem
%%Instance
meanPCG_J_time=zeros(numofProblems,1);
meanPCG_JK_time=zeros(numofProblems,1);
meanPCG_ratio_time=zeros(numofProblems,1);

%%%Save Densities of A matrix
densityA=zeros(numofProblems,1);



%%Needed for sprandsym for generating A matrix
kind=1;


%%%Tolerance and maxiter for PCG
tolPCG = paramsPCG.tolPCG;
maxiterPCG = paramsPCG.maxiterPCG;

trialNum=0;  % trial number

for i=1:numofProblems
    n=dimvec(i);
    omegaa = @(B)( (trace(B)/n)/(det_rootn(B)) );
    seed=seedvec(i);
    rng(seed);
    b=randn(n,1);
    density=densityvec(i);
    if density <= 3/n
        density = density + (2/n);
    end
    rc=rcvec(i);
    %%Create Sparse Positive Definite Matrix using sprandsym
    startgen = tic;
    A=sprandsym(n,density,rc,kind);
    endgen = toc(startgen);

    densityA(i)=nnz(A)/numel(A);
    fprintf('NEW: gentime %g; size n %i; density = %g; condest %g \n', ...
        endgen,n,densityA(i),condest(A))  % use lam1/lamn instead???
    %%%Jacobi Preconditioner
    d=diag(A);
    M=diag(d);
    DJ = diag(sqrt(d));
    %end
    J = DJ\A/DJ;  % Jacobi preconditioning
    J = (J+J')/2;

    %dbJ = DJ\b;

    %%%smallest eigenvalue of J
    lamn_J = eigs(J, 1, 'SM');

    %%%largest eigenvalue of J
    lam1_J= eigs(J, 1, 'LM');


    %%%Gradient descent with respect to kappa on J: MaxIter: %%%100
    startGK = tic;
    DJK=GDKAPPA(paramsGD,J);
    endGK = toc(startGK);
    
    JK=DJK*J*DJK;
    JK=(JK+JK')/2;

    % %%%Testing
    % dJtest=diag(JK);
    % Mtest=diag(dJtest);
    % Dtest = diag(sqrt(dJtest));
    % Jtest = Dtest\JK/Dtest; 
    % 
    % ans1=inv(Dtest)*DJK;
    % ans2=normest(sqrt(diag(JK))-diag(DJK))
    % ans3=normest(ans1-speye(n))
    % ans4=normest(J-Jtest)
    % keyboard

    HK=DJ\DJK;  % change to vectors   dJK./dJ ????
    
    %keyboard 
    %dk=1./(diag(DJK.^2));
    dk=1./(diag(HK.^2));   % use square from start?
    MK=spdiags(dk,0,n,n);

    
    %keyboard
    %MK=DJK^2;
    %dbJK = sqrt(DJK)*b;

    %%%smallest eigenvalue of JK
    lamn_JK = eigs(JK, 1, 'SM');

    %%%Largest eigenvalue of JK
    lam1_JK= eigs(JK, 1, 'LM');


    kappaJ = lam1_J/lamn_J;
    kappaJK = lam1_JK/lamn_JK;

    %condDiff(probNum)=(lam1_J/lamn_J)-(lam1_JK/lamn_JK);
    condDiff(i)= 100*(kappaJK-kappaJ)/(1+((kappaJK-kappaJ)/2));
    condRatios(i)= kappaJK/kappaJ;

    omegaaJK = omegaa(JK);
    omegaaJ = omegaa(J);

    %omegaDiff(probNum)=omegaa(J)-omegaa(JK);
    omegaDiff(i) = 100*(omegaaJK-omegaaJ)/(1+((omegaaJ+omegaaJ)/2));
    omegaRatios(i)=omegaaJK/omegaaJ;



    for k=1:numofInitialPoints
        trialNum=trialNum+1;

        rs = RandomStartPointSet(NumStartPoints=1);
        %%%Create Initial Points for PCG
        problem = createOptimProblem("fmincon",x0=zeros(n,1),lb=-Inf*ones(n,1),ub=Inf*ones(n,1));
        x0 = transpose(list(rs,problem));
        %x0=randn(n,1);

        %[xJ,flagJ,relresJ,iterJ]= ...
        %pcg(J,dbJ,tolPCG,maxiterPCG,[],[],x0);

        startJpcg = tic;
        [xJ,flagJ,relresJ,iterJ]=pcg(A,b,tolPCG,maxiterPCG,M,[],x0);
        timeJ=toc(startJpcg);
	if flagJ
		fprintf('flagJ pcg is %i \n',flagJ)
	end

        %%%For debugging
        % test1=inv(sqrt(M))*A*inv(sqrt(M));
        % lamn_check1 = eigs(test1, 1, 'SM');
        % lam1_check1 = eigs(test1, 1, 'LM');
        % kappacheck1=lam1_check1/lamn_check1;

        % timeJ
        % keyboard
        
        % iterJ
        % 
        % keyboard


        %[xJK,flagJK,relresJK,iterJK]=
        %pcg(JK,dbJK,tolPCG,maxiterPCG,[],[],x0);

        startJKpcg = tic;
        [xJK,flagJK,relresJK,iterJK]= ...
            pcg(A,b,tolPCG,maxiterPCG,MK,[],x0);
        timeJK=toc(startJKpcg) + endGK;
	if flagJK
		fprintf('flagJK pcg is %i \n',flagJK)
	end
        
        %keyboard

        %%%% For debugging
        % test2=inv(sqrt(MK))*A*inv(sqrt(MK));
        % lamn_check2 = eigs(test2, 1, 'SM');
        % lam1_check2 = eigs(test2, 1, 'LM');
        % kappacheck2=lam1_check2/lamn_check2;

        %keyboard
     

        % timeJK
        % keyboard

        % iterJK
        % keyboard


        %%%Fill in PCG iterations for each trial before and after
        %%%preconditioning
        PCG_J_iter(trialNum)=iterJ;
        PCG_JK_iter(trialNum)=iterJK;
        PCG_ratio_iter(trialNum)=iterJK/iterJ;

        %%%%Fill in PCG time for each trial before and after
        %%%%preconditioning
        PCG_J_time(trialNum)=timeJ;
        PCG_JK_time(trialNum)=timeJK;
        PCG_ratio_time(trialNum)=timeJK/timeJ;

    end
end
%keyboard

%%%%Compute Mean Number of PCG Iterations and Time For Each Problem Instance Using
%%%%Many Different Starting Points
for i=1:numofProblems
    meanPCG_J_iter(i)=mean(PCG_J_iter((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_JK_iter(i)=mean(PCG_JK_iter((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_ratio_iter(i)=mean(PCG_ratio_iter((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));

    meanPCG_J_time(i)=mean(PCG_J_time((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_JK_time(i)=mean(PCG_JK_time((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_ratio_time(i)=mean(PCG_ratio_time((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));

end

%keyboard

%keyboard

%%%%Create LaTeX Table

fmt1 = '%5.0f & %5.1e &';


fmt3 = '%6.2e & %7.4e & %4.1f & %7.3e & %4.2f & %4.2f';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-2}\\cline{2-8}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|cc|cc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{2}{|c||}{Dim \& Density} & ', ...
    '\multicolumn{2}{|c|}{Ratios Conds $(K,J)$ \#} & ', ...
    '\multicolumn{2}{|c|}{K: $\kappa$-Desc.} &', ...
    '\multicolumn{2}{|c|}{Ratios $K/J$}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density');
fprintf(fid, hfmt3,'$\kappa(K)/\kappa(J)  $', '$\omega(K)/\omega(J)  $', ...
    'iters', 'cpu', ...
    'iters', 'cpu');
fprintf(fid, fmt4);

for k=1:numofProblems
    fprintf(fid, fmt1, dimvec(k), densityA(k));
    fprintf(fid, fmt3, condRatios(k), omegaRatios(k), ...
	  meanPCG_JK_iter(k), meanPCG_JK_time(k), ...
            meanPCG_ratio_iter(k), meanPCG_ratio_time(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);




end













%%%Start Gradient Descent Function. Takes in Parameters PAR and Input
%%%Matrix/Initial Matrix A;
function D_update_product=GDKAPPA(PAR,Ainput)
tic
n=size(Ainput,2);
identity=speye(n);
[xn, lamn] = eigs(Ainput, 1, 'SM');
%%%Largest eigenvalue and eigenvector
[x1, lam1] = eigs(Ainput, 1, 'LM');

%%%% Create Kappa descent direction from eigenvectors of A.
gradient=(((lamn*lam1)*(x1.*x1))-((lam1*lamn)*(xn.*xn)));
stopcrit=norm(gradient);
gradient(abs(gradient)<1e-16)=0;
gradientMatrixsparse=spdiags(gradient,0,n,n);

%%Iteration count and flag
iter=1;
outeriter=1;
%D_update=identity;
D_update_product=identity;

if stopcrit<=PAR.tol
    A_update=Ainput;
    D_update_product=identity;
    KappaA_update=lam1/lamn;
end

%%Initialize A which is continuously updated in algorithm
A=Ainput;
while outeriter<=PAR.maxiter
    %while stopcrit>PAR.tol
    Kappa_A=lam1/lamn;

    D_update=identity-((1/PAR.L)*gradientMatrixsparse);

    sqrtD_update=sqrt(D_update);


    % norm(diag(D_update))
    A_update=sqrtD_update*A*sqrtD_update;
    A_update=(A_update+A_update')/2;

    D_update_product_old=D_update_product;
    D_update_product=(sqrtD_update)*D_update_product_old;
   

    % keyboard
    % B=D_update_product*Ainput*D_update_product;
    % normdiff=normest(B-A_update)
    % keyboard



    [xn, lamn, flag1] = eigs(A_update, 1, 'SM');
    %%%Largest eigenvalue and eigenvector
    [x1, lam1, flag2] = eigs(A_update, 1, 'LM');

    
    KappaA_update=lam1/lamn;


    NormDiff=(PAR.norm_fn(D_update-identity))^2;
    Ftest = Kappa_A+PAR.prod_fn(D_update-identity,gradientMatrixsparse)+(0.5*PAR.L*NormDiff)+10^-6-KappaA_update;

    if flag1==1 || flag2==1
        A_update=A;
        D_update_product=D_update_product_old;
        break
    end



    %%%%Backtracking Line Search for Lipschitz
    while Ftest < 0
        PAR.L=PAR.beta*PAR.L;
        %%%%% update
        Kappa_A=lam1/lamn;

        %%%Gradient Descent Step
        D_update=identity-((1/PAR.L)*gradientMatrixsparse);
        A_update=sqrt(D_update)*A*sqrt(D_update);
        A_update=(A_update+A_update')/2;

        [xn, lamn] = eigs(A_update, 1, 'SM');
        %%%Largest eigenvalue and eigenvector
        [x1, lam1] = eigs(A_update, 1, 'LM');

        KappaA_update=lam1/lamn;
        NormDiff=(PAR.norm_fn(D_update-identity))^2;

        Ftest = Kappa_A+PAR.prod_fn(D_update-identity,gradientMatrixsparse)+(0.5*PAR.L*NormDiff)+10^-6-KappaA_update;

        %%Compute Gradient
        gradient=(((lamn*lam1)*(x1.*x1))-((lam1*lamn)*(xn.*xn)));
        iter=iter+1;
        stopcrit=norm(gradient);
        PAR.L;
        gradient(abs(gradient)<1e-16)=0;
        gradientMatrixsparse=spdiags(gradient,0,n,n);



        A=A_update;



    end
    %%%Compute Gradient
    gradient=(((lamn*lam1)*(x1.*x1))-((lam1*lamn)*(xn.*xn)));
    iter=iter+1;
    outeriter=outeriter+1;
    %%Stopcrit is Norm of Gradient
    stopcrit=norm(gradient);
    gradient(abs(gradient)<1e-16)=0;
    gradientMatrixsparse=spdiags(gradient,0,n,n);
    A=A_update;

    % if mod(iter,20) == 0
    %     fprintf('total iter = %i; Normofgrad = %g; Condition Number Value=%g\n', ...
    %         iter, stopcrit,  KappaA_update)
    % end


end
totalTime=toc;
end

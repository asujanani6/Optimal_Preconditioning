function run_testsOptKappa_Sparse_PCG(dimvec, seedvec, densityvec, numofInitialPoints, params, filename)
%%% Start with a kappa-optimal scaled A and compare to apply 
%%% omega scaling to it

%%%Uses function  that is below
%% function [A,b]=findoptkappaASparse(n,seed,density)


%%    det_rootn(B) is used below so cvx is needed
if ~exist('det_rootn')
	cvx_setup
end

%%%Number of Problems is just number of dimensions inputted. Number of
%%%dimensions should be same as number of seeds and number of densities
numofProblems=length(dimvec);

%%Number of Trials (with Many Different Initial Points) and Number of
%%Problems
numofTrials=numofInitialPoints*numofProblems;


%%%%% vectors to save output. The difference
%%%%% between classical condition number of optimal kappa and Jacobi is also saved
condDiff=zeros(numofProblems,1);
condRatios=zeros(numofProblems,1);

%For each problem, the difference
%%%%% between omega condition number of original matrix and Jacobi is saved
omegaDiff=zeros(numofProblems,1);
omegaRatios=zeros(numofProblems,1);

%%%Save PCG iterations before and after Preconditioning as Well as Ratios
PCG_DAD_iter=zeros(numofTrials,1);
PCG_A_iter=zeros(numofTrials,1);
PCG_ratio_iter=zeros(numofTrials,1);


%%%Save PCG time before and after Preconditioning as Well as Ratios
PCG_DAD_time=zeros(numofTrials,1);
PCG_A_time=zeros(numofTrials,1);
PCG_ratio_time=zeros(numofTrials,1);


%%Save Mean Iterations over many different Initial Points for Each Problem
%%Instance
meanPCG_DAD_iter=zeros(numofProblems,1);
meanPCG_A_iter=zeros(numofProblems,1);
meanPCG_ratio_iter=zeros(numofProblems,1);

%%Save Mean Time over many different Initial Points for Each Problem
%%Instance
meanPCG_DAD_time=zeros(numofProblems,1);
meanPCG_A_time=zeros(numofProblems,1);
meanPCG_ratio_time=zeros(numofProblems,1);

%%%Save Densities of A matrix
densityA=zeros(numofProblems,1);

% %%%Save dimensions of A matrix
% % dimAplace=zeros(numofTrials,1);
% dimA=zeros(numofProblems,1);


tol = params.tol;
maxiter = params.maxiter;


trialNum=0;  % trial number


for i=1:numofProblems
    n=dimvec(i);
    omegaa = @(B)( (trace(B)/n)/(det_rootn(B)) );
    seed=seedvec(i);
    density=densityvec(i);
    startgen = tic;
    [A,b]= findoptkappaASparse(n,seed,density);
    endgen = toc(startgen);
    %%%Save densities of A as a vector
    densityA(i)=nnz(A)/numel(A);
    fprintf('NEW: time gen = %g; size n =%i; density %g density-2/n %g \n', ...
        endgen,n,densityA(i),densityA(i)-(2/n))
    %%%Jacobi Preconditioner
    %D = diag(sqrt(diag(A)));
    d = diag(A);
    M = diag(d);
    D = diag(sqrt(d));
    %end
    DAD = D\A/D;  % Jacobi preconditioning
    DAD = (DAD+DAD')/2;
    %db = D\b;

    %%%smallest eigenvalue of DAD
    lamn_DAD = eigs(DAD, 1, 'SM');

    %%%Largest eigenvalue of DAD
    lam1_DAD= eigs(DAD, 1, 'LM');

    %%%smallest eigenvalue of A
    lamn_A = eigs(A, 1, 'SM');

    %%%Largest eigenvalue of A
    lam1_A= eigs(A, 1, 'LM');

    % condDiff(probNum)=(lam1_A/lamn_A)-(lam1_DAD/lamn_DAD);
    % omegaDiff(probNum)=omegaa(A)-omegaa(DAD);

    condDAD = lam1_DAD/lamn_DAD;
    condA = lam1_A/lamn_A;
    condDiff(i)= ...
        100*(condDAD-condA)/(1+((condDAD+condA)/2));
    condRatios(i)= ...
        condDAD/condA;

    omegaaDAD=omegaa(DAD);
    omegaaA=omegaa(A);
    omegaDiff(i)=100*(omegaaDAD-omegaaA) ...
        /(1+((omegaaDAD+omegaaA)/2));  % for percent change in J
    omegaRatios(i)=omegaaDAD/omegaaA;

    for k=1:numofInitialPoints
        trialNum=trialNum+1;
        rs = RandomStartPointSet(NumStartPoints=1);
        %%%Create Initial Points for PCG
        problem = createOptimProblem("fmincon", ...
            x0=zeros(n,1),lb=-Inf*ones(n,1),ub=Inf*ones(n,1));
        x0 = transpose(list(rs,problem));
        %x0=randn(n,1);

        %[xA,flagA,relresA,iterA]= pcg(A,b,tol,maxiter,[],[],x0);

        tic
        [~,flagA,~,iterA]= pcg(A,b,tol,maxiter,[],[],x0);
        timeA=toc;
        if flagA ~= 0
            fprintf('pcg A failed\n')
            % keyboard
        end

        %[xDAD,flagDAD,relresDAD,iterDAD]= ...
        %pcg(DAD,db,tol,maxiter,[],[],x0);

        tic
        [~,flagDAD,relresDAD,iterDAD]= ...
            pcg(A,b,tol,maxiter,M,[],x0);  %D=sqrt M
        timeDAD=toc;

        % if flagDAD ~= 0
        % 	fprintf('DAD failed; flag,relres,iter= \n', ...
        %         flagDAD,relresDAD,iterDAD)
        % 	keyboard
        % end

        %%%Fill in PCG iterations for each trial before and after
        %%%preconditioning
        PCG_DAD_iter(trialNum)=iterDAD;
        PCG_A_iter(trialNum)=iterA;
        PCG_ratio_iter(trialNum)=iterA/iterDAD;

        %%%%Fill in PCG time for each trial before and after
        %%%%preconditioning
        PCG_DAD_time(trialNum)=timeDAD;
        PCG_A_time(trialNum)=timeA;
        PCG_ratio_time(trialNum)=timeA/timeDAD;

    end
end

%%%%Compute Mean Number of PCG Iterations and Time For Each Problem Instance Using
%%%%Many Different Starting Points
for i=1:numofProblems
    meanPCG_DAD_iter(i)=mean(PCG_DAD_iter((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_A_iter(i)=mean(PCG_A_iter((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_ratio_iter(i)=mean(PCG_ratio_iter((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));

    meanPCG_DAD_time(i)=mean(PCG_DAD_time((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_A_time(i)=mean(PCG_A_time((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));
    meanPCG_ratio_time(i)=mean(PCG_ratio_time((numofInitialPoints*(i-1)+1):(numofInitialPoints*i)));


end


%keyboard


%%%%Create LaTeX Table

fmt1 = '%6.0f & %6.2e &';


fmt3 = '%6.2e & %7.3e & %4.1f & %7.3e & %3.1f & %3.1f';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-2}\\cline{2-8}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|cc|cc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{2}{|c||}{Dim \& Density} & ', ...
    '\multicolumn{2}{|c|}{Ratios in conds $(J,A)$ } & ', ...
    '\multicolumn{2}{|c|}{$J$: $\omega$-opt of $A$} &', ...
    '\multicolumn{2}{|c|}{Ratios A/J}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$n$', 'density');
fprintf(fid, hfmt3,'$\kappa(J)/\kappa(A)  $', '$\omega(J)/\omega(A)  $', ...
    'iters', 'cpu', ...
    'iters', 'cpu');
fprintf(fid, fmt4);

       %meanPCG_DAD_iter(k), meanPCG_A_iter(k), ...
for k=1:numofProblems
    fprintf(fid, fmt1, dimvec(k), densityA(k));
    fprintf(fid, fmt3, condRatios(k), omegaRatios(k), ...
       meanPCG_DAD_iter(k), meanPCG_DAD_time(k), ...
           meanPCG_ratio_iter(k), meanPCG_ratio_time(k));
    fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);




end



%%%Sparse Version for Finding Optimal Kappa A
function [A,b]=findoptkappaASparse(n,seed,density)
rng(seed)
if exist('n','var')
    if mod(n,2) ~=0
        n = 2*round(n/2);
    end
else
    n = 1000;
end
%density0=0.001
%fprintf('\nNEW PROB n  %i \n',n)


lam1 = 1e2*(rand+.5);   % large   ?????eventually this is input???
lamn = 1e-2*(rand+.5);   % small   ?????eventually this is input???
k = round(n/2);
%%%%%%%%%generate x1,xn orthog and hadmard equal
x0 = sprandn(k,1,density);  % use half the vector with different sign
x0 = x0/norm(x0);
x00 = sprandn(k,1,density);  % use half the vector with different sign
x00 = x00/norm(x00);
xn = [x0; -x00];
xn = xn/norm(xn);
x1 = [x0; x00];
x1 = x1/norm(x1);
rperm = randperm(n);
x1 = x1(rperm);
xn = xn(rperm);

% x1'*xn;
% norm(x1.*x1-xn.*xn);
%%% fix check below????????
%fprintf('check orthog x1''xn and hadamard diff orthog %e %e\n', ...
%         x1'*xn,norm(x1.*x1-xn.*xn))  % a check on had. orthog.
%Q = null([x1 xn]');   % remaining n-2 eigenvectors

X1n = [x1 xn];
[qq,~] = qr(X1n);   % get orthon basis for [x1 xn]^perp
Q = qq(:,3:end);
norm(X1n'*Q,'fro');
meanlam = mean([lam1 lamn]);
%lam = meanlam + ((meanlam/4)*((rand(n-2,1)+.5)-1) );
lam = linspace(lamn+0.0001,lam1-3,n-2);

%%%Use spdiags so that Amid is in sparse format
Amid = Q*spdiags(lam',0,n-2,n-2)*Q'; %%Construct A mid using middle eigenvalues.

%%%Remove very small numbers less than machine precision
%Amid(abs(Amid)<1e-15)=0;  % this takes 90% of the time????????
%Amid = (Amid+Amid')/2; %symmetrize to avoid roundoff

A = lam1*(x1*x1')+lamn*(xn*xn')+Amid;
A = (A+A')/2;    % A is optimal

b = randn(n,1);
end

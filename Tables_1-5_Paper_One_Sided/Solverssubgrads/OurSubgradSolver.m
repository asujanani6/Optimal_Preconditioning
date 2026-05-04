function output = OurSubgradSolver(A,params,options)
%function output = HenryOptKappa(A,params,options)
% Input:   A n by n pos def
%          params and options see below for defaults
%% min kappa(A(e+Vv)), where cols of V is basis for e^perp.
%%%  using eigs for min/max eigs and ignoring possible multiplicities
%%% therefore the algorithm is ignoring the nondiff/nonsmmooth of kappa
%%%comments: Aug7 updates polished up to NOT repeat unecessarily
%%%Aug9: need to replace eigs by e.g. max eig  A(Diag(w)) and test
%%% the roundoff error problem???

%%% params structure with defaults
if isfield(params,'R')  % Cholesky for A
	R = params.R;
else
	R = chol(A);  % A = R'*R
end
if isfield(params,'maxitermain')  % iter bound for main while loop
	maxitermain = params.maxitermain;
else
	maxitermain = 50;
end

n=length(A);
if isfield(params,'steplgthmin')  % stop bnd for steps are too small
	steplgthmin = params.steplgthmin;
else
	steplgthmin = 1e-7/n;
end

if isfield(params,'tolerance')  % tolerance for main while loop
	tolerance = params.tolerance;
else
	tolerance = 1e-3;
end


if isfield(params,'x1')  % for warm starts
	x1 = params.x1;
else
	x1 = randn(n,1);
	x1 = x1/norm(x1);
end


if isfield(params,'xn')  % for warm starts
	xn = params.xn;
else
	xn = randn(n,1);
	xn = xn/norm(xn);
end

%%% options structure with defaults
if isfield(options,'plots')  % plots for debugging mainly
	plots = options.plots;
else
	plots = false;
end



%%%Generate starting points
sqrt2 = sqrt(2); % for fast V v mult.
%e = ones(n,1);
%V = sparse([speye(n-1); -ones(1,n-1)]/sqrt2);  % changed to efficient use
v=zeros(n-1,1);
w = 1 + [v;-sum(v)]/sqrt2;   %1 + V*v;
%%%Maximum Eigenpair for A*Diag(e+V*v)  *x  = A* ( (e+V*v).*x)
%[R,flag,P] = chol(S) additionally returns a permutation matrix P, which is a preordering of sparse matrix S obtained by amd. If flag = 0, then S is symmetric positive definite and R is an upper triangular matrix satisfying R'*R = P'*S*P.
%%%%% iteration 1 before entering MAIN while loop
Afunctionalmax = @(x)( A*(w.*x) );   % for max eig  A(Diag(w))
% Afunctionalmin = @(x)( (R\(R'\(x./w))) );   % for min eig
Afunctionalmin = @(x) (R\(R'\x))./w;

[x1, lam1,flag1]= eigs(Afunctionalmax ,n,1, 'LM', ...
	       'StartVector', x1, 'Tolerance', 1e-14);
[xn, lamn,flagn] = eigs(Afunctionalmin, n, 1, 'SM', ...
	       'StartVector', xn, 'Tolerance', 1e-14);
kappa = lam1/lamn;

if flag1==1 || flagn==1
 fprintf("eigenvector computation failed")
end
%%%%%%%%%%%%%
% gtemp = ((x1.*x1)/norm(w.*x1)-(xn.*xn)/norm(w.*xn)); 
gtemp = ((x1.*x1)/dot(x1,w.*x1)-(xn.*xn)/dot(xn,w.*xn));
gradv = kappa*(gtemp(1:n-1) -gtemp(n))/sqrt2;  % kappa*V'*.....
%gradvold = V'*((x1.*x1)/norm(w.*x1)-(xn.*xn)/norm(w.*xn)); % without kappa
deltav = -gradv; % direction of change; steepest descent
tstep = inf; % no first steplength


relnormg=Inf;  % for stopping main while loop
noiter = 1;

%%%MAIN WHILE
while relnormg>tolerance && noiter < maxitermain  && tstep > steplgthmin
	

	%%% save for output
	relnormgs(noiter) = relnormg;
	kappas(noiter) = kappa;
	tsteplgths(noiter) = tstep;
	%%%%%%%%%%%%%
	
	%% prepare for linesearch for steplengh e + Vv + tV deltav > 0
	%%   w = e+Vv;  w+tV deltav > 0;  w + t deltaw > sigma w, 0<sigma small
	%% sigma = .03; wi + t*deltawi > sigma wi; t*deltawi > (sigma-1)wi 
	%% t*deltawi > -.97wi, deltawi < 0;  wbar = -.97w
	%% deltawi<0: t < wbari./deltawi   deltawi < 0 
	%%% .9 to stay away from boundary; should .97 be better???
	wbar = .97*w;  % wbar > 0  and uses the latest v
	%deltaw = V*deltav;
	deltaw = [deltav;-sum(deltav)]/sqrt2;
	%tmax = min( 1, ...
	%	min((-wbar(deltaw<0))./deltaw(deltaw<0)) ); %safeguard to bdry
	tmax = min((-wbar(deltaw<0))./deltaw(deltaw<0)) ; %safeguard to bdry

	%%%%%%%%%%%%%%%%%for debugging%%%%%%%%%%%%%%%%%%
	if plots  % plot kappa and directderiv
		ts = linspace(-tmax/10,0,8);
		ts = [ts linspace(1e-3*tmax,tmax,8)];
		kappats = [];
		dirderts = [];
		lamsvecs = [];
		itertt = 0;
		multeigs = 1;
		for tt = ts
			itertt = itertt + 1;
			vtest=v+tt*deltav;  % cheap linesearch and take step
			%wts=e+V*vtest;
			wts=1+([vtest;-sum(vtest)]/sqrt2);
			ADw = A*diag(wts);
			[X, Lam] = eig(ADw);
			[lams,indslam] = sort(diag(Lam),'descend');
			lamsvecs = [lamsvecs lams];
			lamn = sum(lams(end-multeigs+1:end))/multeigs; 
			lam1 = lams(1)/multeigs; 
			X1 = X(:,indslam(1));
			Xn = X(:,indslam(end));
			%lam1 = sum(lams(1:multeigs))/multeigs; 
			kappats(itertt) = lam1/lamn;
			%gtemp = ((X1.*X1)/norm(wts.*X1)-(Xn.*Xn)/norm(wts.*Xn));
            gtemp = ((x1.*x1)/dot(X1,wts.*X1)-(Xn.*Xn)/dot(Xn,wts.*Xn));
			gradvtest = ...
			    kappats(itertt)*(gtemp(1:n-1) -gtemp(n))/sqrt2;
			%gradvtest = kappats(itertt)*V'*((X1.*X1)/ ...
			%	norm(wts.*X1)-(Xn.*Xn)/norm(wts.*Xn));
			dirderts(itertt) = deltav'*gradvtest;
			%[kappats(1:itertt)' dirderts(1:itertt)']
			%figure(2)
			%clf
			%semilogy(diff(lamsvecs))
			%drawnow
	 	end
			figure(2)
			clf
			plot(ts(1:itertt),kappats,'-x')
			title(['ts versus kappa; iter = ',num2str(noiter)])
			drawnow
			figure(3)
			clf
			plot(ts(1:itertt),dirderts,'-o')
			title(['ts versus dirderiv; iter = ',num2str(noiter)])
			drawnow
	%%??????temp test????????????????????????
	end
	%%%%%%%%%%%%%%%%%end for debugging%%%%%%%%%%%%%%%%%%


	% linesearch using directional derivative/backtracking
	%% prepare for first backtrack before entering while backtrack loop
	vtest=v+tmax*deltav;
	%%% need wbar + t deltaw > 0;  t deltaw > -wbar
	%%%    t <=  -wbar_i/deltaw_i    deltaw_i < 0
	%wtest=e+V*vtest;
	wtest=1+([vtest;-sum(vtest)]/sqrt2);
	Afunctionalmax = @(x)( A*(wtest.*x) );   % for max eig
	%Afunctionalmin = @(x)( (R\(R'\(x./wtest))) );   % for min eig
    Afunctionalmin = @(x) (R\(R'\x))./wtest;

	[x1, lam1, flag1]= eigs(Afunctionalmax ,n,1, 'LM', ...
		       'StartVector', x1, 'Tolerance', 1e-14);
    	[xn, lamn, flagn] = eigs(Afunctionalmin, n, 1, 'SM', ...
		       'StartVector', xn, 'Tolerance', 1e-14);
	kappatest = lam1/lamn;

    if flag1==1 || flagn==1
        fprintf("eigenvector computation failed")
    end

	% gtemp = ((x1.*x1)/norm(w.*x1)-(xn.*xn)/norm(w.*xn)); 
    gtemp = ((x1.*x1)/dot(x1,w.*x1)-(xn.*xn)/dot(xn,w.*xn));

	gradvtest = kappatest*(gtemp(1:n-1) -gtemp(n))/sqrt2;
	%gradvtest = kappatest*V'*((x1.*x1) ...
	%	/norm(wtest.*x1)-(xn.*xn)/norm(wtest.*xn));
	slopetest = deltav'*gradvtest;  % directional deriv in dir deltav

	iterlnsrch = 1;
	tstep = tmax;  % starting steplength guarantees pos def
	while (slopetest > 0 || kappatest > kappa)   ...
		     && tstep > steplgthmin 
		tstep = tstep/2;   % backtrack 'a lot'
		vtest=v+tstep*deltav;  % backtrack
		%%% need wbar + t deltaw > 0;  t deltaw > -wbar
		%%%    t <=  -wbar_i/deltaw_i    deltaw_i < 0
		%wtest=e+V*vtest;
		wtest=1+([vtest;-sum(vtest)]/sqrt2);
		Afunctionalmax = @(x)( A*(wtest.*x) );   % for max eig
		%Afunctionalmin = @(x)( (R\(R'\(x./wtest))) );   % for max eig
        Afunctionalmin = @(x) (R\(R'\x))./wtest;

		%[x1, lam1]= eigs(Afunctionalmax ,n,1);
		%[xn, lamn] = eigs(Afunctionalmin, n, 1, 'SM');
		[x1, lam1, flag1]= eigs(Afunctionalmax ,n,1, 'LM', ...
			       'StartVector', x1, 'Tolerance', 1e-14);
	    	[xn, lamn, flagn] = eigs(Afunctionalmin, n, 1, 'SM', ...
			       'StartVector', xn, 'Tolerance', 1e-14);
		kappatest = lam1/lamn;

        if flag1==1 || flagn==1
            fprintf("eigenvector computation failed")
        end
	
		%gtemp = ((x1.*x1)/norm(w.*x1)-(xn.*xn)/norm(w.*xn)); 
        gtemp = ((x1.*x1)/dot(x1,w.*x1)-(xn.*xn)/dot(xn,w.*xn));
		gradvtest = kappatest*(gtemp(1:n-1) -gtemp(n))/sqrt2; %kappa*V'*...
		%gradvtest = kappatest*V'*((x1.*x1)/ ...
		%	 norm(wtest.*x1)-(xn.*xn)/norm(wtest.*xn));
		%deltavtest = -gradvtest;  % direction steepest descent
		slopetest = deltav'*gradvtest; % slope at new point on line
		iterlnsrch = 1 + iterlnsrch;

	end


	% UPDATES if successful
	if  tstep <= steplgthmin % line search failed
		tstep = 0;
		%kappatest = kappa;  % unchanged
		%gradvtest = gradv; % unchanged
	else
		v = vtest; % successful linesearch with decrease
		kappa = kappatest;
		gradv = gradvtest;
		deltav = -gradv;  % direction of change; steepest descent
		w=1+([v;-sum(v)]/sqrt2);
		relnormg = (norm(gradv)/kappa)^2; % relnorm for line search
	end



	%%% need wbar + t deltaw > 0;  t deltaw > -wbar
	%%%    t <=  -wbar_i/deltaw_i    deltaw_i < 0
	%w=e+V*v;
	
	
	noiter = noiter + 1;
	relnormgs(noiter) = relnormg;
	kappas(noiter) = kappa;

	if plots
		figure(1)
		clf
		semilogy(kappas,'-x')
		%hold on
		%semilogy(relnormgs)
		%legend('kappa','relnorm','location','best')
		title('iter vs kappa')
		drawnow
	
	

	end  % of if plots
	
	

end  % of main while



output.v = v;
output.noiter = noiter;
output.kappas = kappas;
output.relnormgs= relnormgs;
output.tsteplgths= tsteplgths;
output.w=w;


end  % of function







function [csqs,dsqs,output] = Solver_Alt_TwoSided(A,optionsOurSolver)

verbose=optionsOurSolver.verbose;

%% we have to assume NO zero cols/rows
%%% check to return if a zero row or col
indscols = find(all(A==0));
indsrows = find(all(A == 0, 2));
if length(indscols)+length(indsrows) > 0
	fprintf('ERROR IN DATA: zero rows or cols; rerun for new problem\n')
	return  
end


dsqone=zeros(optionsOurSolver.maxiterscale,1);
csqone=zeros(optionsOurSolver.maxiterscale,1);
dsqone(1)=inf;
csqone(1)=inf;
n=size(A,2);

csqs = ones(n,1);   % for multiplications
dsqs = ones(n,1);   % for multiplications
omegaas = 0;
iter = 0;
maxnormones = inf;
%%%% alternate col/row scalings
tic
while maxnormones > optionsOurSolver.tolcd && iter <= optionsOurSolver.maxiterscale
	iter = iter + 1;
	%%% vector cnvges to all ones so norm will be sqrt n
	dsqinv = vecnorm(A)'; %  inv sqrt Diag(d) in paper
    % dsqinv=sqrt(vecnorm(A)');
	dsqone(iter+1) = norm(dsqinv-1);
	A = A*spdiags(1./dsqinv,0,n,n);
	dsqs = dsqs.*dsqinv;  % for diag(M2)
	%%% vector cnvges to all ones so norm will be sqrt n
	csqinv = vecnorm(A,2,2); %  inv sqrt Diag(c) in paper
    %csqinv = sqrt(vecnorm(A,2,2));
    csqone(iter+1)=norm(csqinv-1);
	A = spdiags(1./csqinv,0,n,n)*A;
	csqs = csqs.*csqinv;  % for diag(M1)
	if verbose
        	omegaas(iter) = omegaa(A);
	else
        	omegaas(iter) = 0;
	end
	maxnormones = ...
            max(abs(diff(dsqone(iter:iter+1))), ...
                    abs(diff(csqone(iter:iter+1))));  % max changes
    % keyboard
end
output.time=toc;
output.iter=iter;




fprintf('\n\nNEW: while ends with: max(norm(c-1),norm(d-1)) %g  iter %i\n', ...
         maxnormones,iter)


if verbose
	Mcd = @(c,d)( A'*spdiags(c,0,n,n)*A*spdiags(d,0,n,n) );
	fcd = @(c,d)( trace(Mcd(c,d))/n );
	gcd = @(c,d)( det_rootn(Mcd(c,d)) );
	Gcd = @(c,d)([ diag(A*spdiags(d,0,n,n)*A') 
               diag(A'*spdiags(c,0,n,n)*A)] );
	Fcd = @(c,d)(Gcd(c,d) - fcd(c,d)*[1./c; 1./d ] );
	Fcdcd = Fcd(csqs,dsqs);
	fprintf('norm F(c,d) G-fe for M1,M2 =%g\n',norm(Fcdcd,'fro'));
end

if verbose
	%% plot from point 2 on  as first point has inf values
	figure(1)
	clf
	plot(csqone(2:end))
	hold on
	plot(dsqone(2:end))
	plot(omegaas(1:end))
	legend('csqone','dsqone','omegaas','location','best')
	hold off
	drawnow
end


end






   
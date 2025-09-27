function output = OurSubgradSolverSimplex(A,params)
% Input:   A n by n pos def
%          params and options see below for defaults
%% min kappa(A(e+Vv)), where cols of V is basis for e^perp.

%%% params structure with defaults
if isfield(params,'R')  % Cholesky for A
    R = params.R;
else
    R = chol(A);  % A = R'*R
end
if isfield(params,'maxitermain')  % iter bound for main while loop
    maxitermain = params.maxitermain;
else
    maxitermain = 500;
end

if isfield(params,'hatdelta')  % hatdelta
    hatdelta = params.hatdelta;
else
    hatdelta=1e-3;
end

n=length(A);

if isfield(params,'tolerance')  % tolerance for main while loop
    tolerance = params.tolerance;
else
    tolerance = 1e-4;
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
%Afunctionalmin = @(x)( (R\(R'\(x./w))) );   % for min eig
Afunctionalmin = @(x) (R\(R'\x))./w;

[x1, lam1, flag1]= eigs(Afunctionalmax ,n,1, 'LM', ...
    'StartVector', x1, 'Tolerance', 1e-14);
[xn, lamn, flagn] = eigs(Afunctionalmin, n, 1, 'SM', ...
    'StartVector', xn, 'Tolerance', 1e-14);
kappa = lam1/lamn;

if flag1==1 || flagn==1
    fprintf("eigenvector computation failed")
end
%%%%%%%%%%%%%

avec=ones(n-1,1);
bval=-sqrt2*(hatdelta-1);
lbS=sqrt2*(hatdelta-1);
ubS=Inf;

noiter = 0;

stopcrit=Inf;

%%%MAIN WHILE
while noiter < maxitermain && stopcrit>tolerance
    noiter = noiter + 1;
    % t=(0.25)/sqrt(noiter);
    t=(1)/sqrt(noiter);

    %gtemp = ((x1.*x1)/norm(w.*x1)-(xn.*xn)/norm(w.*xn));
    gtemp = ((x1.*x1)/dot(x1,w.*x1)-(xn.*xn)/dot(xn,w.*xn));

    gradv = kappa*((gtemp(1:n-1) -gtemp(n))/sqrt2);
    %gradv = kappa*gradv;      % gradient at v
    deltav = gradv/norm(gradv); % direction of change; steepest descent


    v=proj_gen_simp(v-(t*deltav),avec,bval,lbS,ubS);




    w = 1 + [v;-sum(v)]/sqrt2;

    Afunctionalmax = @(x)( A*(w.*x) );   % for max eig
    %Afunctionalmin = @(x)( (R\(R'\(x./w))) );   % for min eig
    Afunctionalmin = @(x) (R\(R'\x))./w;

    [x1, lam1, flag1]= eigs(Afunctionalmax ,n,1, 'LM', ...
        'StartVector', x1, 'Tolerance', 1e-14);

    [xn, lamn, flagn] = eigs(Afunctionalmin, n, 1, 'SM', ...
        'StartVector', xn, 'Tolerance', 1e-14);

    if flag1==1 || flagn==1
        fprintf("eigenvector computation failed")
    end

    %%Kappa Value
    kappaold=kappa;
    kappa = lam1/lamn;

    stopcrit=(abs(kappaold-kappa))/(0.5*(kappaold+kappa));
end





output.v = v;
output.noiter = noiter;
output.w=w;


    function out = proj_gen_simp(x,a,b,lb,ub)
        %proj_gen_simp computes the orthogonal projection onto general
        %simplex
        %                   {x : <a,x> <= b, lb<=x<=ub}
        %
        %  Usage:
        %  out = proj_gen_simp(x,a,b,[lb],[ub])
        %  ===========================================
        %  Input:
        %  x - point to be projected (vector/matrix)
        %  a - vector/matrix
        %  b - scalar
        %  lb - lower bound (vector/matrix/scalar) [default: -inf]
        %  ub - upper bound (vector/matrix/scalar) [default: inf]
        %  ===========================================
        %  Assumptions:
        %  general simplex is nonempty
        %  ===========================================
        %  Output:
        %  out - projection vector

        % This file is part of the FOM package - a collection of first order methods for solving convex optimization problems
        % Copyright (C) 2017 Amir and Nili Beck
        %
        % This program is free software: you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation, either version 3 of the License, or
        % (at your option) any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program.  If not, see <http://www.gnu.org/licenses/>.

        %reading the user x and setting defalut values when required.
        if (nargin < 3)
            error ('usage: proj_halfspace_box(x,a,b,[lb],[ub]') ;
        end

        if (nargin < 5)
            %upper bound is not given, setting to defalut value :inf
            ub = inf ;
        end
        if ((nargin < 4) || (isempty( lb)))
            %lower bound is not given, setting to defalut value :0
            lb = -inf;
        end

        if  (any(any(lb > ub)))
            error('Set is infeasible') ;
        end

        if trace(a'*min(max(x,lb),ub)) <= b
            out =  min(max(x,lb),ub) ;
        else
            %use proj_hyperplane_box
            out = proj_hyperplane_box(x,a,b,lb,ub) ;
        end
    end



    function out = proj_hyperplane_box(x,a,b,l,u)
        %PROJ_HYPERPLANE_BOX computes the orthogonal projection onto the intersection of a hyperplane and a box {x:<a,x>=b,l<=x<=u}
        %
        %  Usage:
        %  out = PROJ_HYPERPLANE_BOX(x,a,b,[l],[u])
        %  ===========================================
        %  Input:
        %  x - point to be projected (vector/matrix)
        %  a - vector/matrix
        %  b - scalar
        %  l - lower bound (vector/matrix/scalar) [default: -inf]
        %  u - upper bound (vector/matrix/scalar) [default: inf]
        %  ===========================================
        %  Assumptions:
        %  The intersection of the hyperplane and the box is nonempty
        %  ===========================================
        %  Output:
        %  out - projection vector

        % This file is part of the FOM package - a collection of first order methods for solving convex optimization problems
        % Copyright (C) 2017 Amir and Nili Beck
        %
        % This program is free software: you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation, either version 3 of the License, or
        % (at your option) any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program.  If not, see <http://www.gnu.org/licenses/>.

        if (nargin < 3)
            error ('usage: proj_hyperplane_box(x,a,b,[l],[u])') ;
        end

        if (nargin < 5)
            %upper bound is not given, setting to defalut value :inf
            u = inf ;
        end
        if ((nargin < 4) || (isempty( l)))
            %lower bound is not given, setting to defalut value :-inf
            l = -inf;
        end

        %checking that sum (l) < b < sum (u)


        sumlb = trace(a'*((l .* (sign(a)>0)) + (u .* (sign(a)<0)))) ;
        sumub = trace(a'*((u .* (sign(a)>0)) + (l .* (sign(a)<0)))) ;

        if ((sumlb > b) || (any(any(l > u))) ||  (sumub < b))
            error('Set is infeasible') ;
        end

        %solve with equality
        %defining f on lambda - to be used by the bisetion

        eps = 1e-10 ;

        f= @(lam)   trace(a'*min(max(x-lam*a,l),u))-b;
        lambda_min = -1;
        while(f(lambda_min)<0)
            lambda_min = lambda_min *2  ;
        end

        lambda_max = 1;
        while(f(lambda_max)>0)
            lambda_max = lambda_max *2  ;
        end

        final_lam = bisection(f,lambda_min,lambda_max,eps) ;
        out=  min(max(x-final_lam*a,l),u) ;
    end


    function out=bisection(f,min_val,max_val,eps)
        %INPUT
        %================
        %f ................... a scalar function
        %lb ................. the initial lower bound
        %ub ................ the initial upper bound
        %eps .............. tolerance parameter
        %OUTPUT
        %================
        % z ................. a root of the equation f(x)=0

        % This file is part of the FOM package - a collection of first order methods for solving convex optimization problems
        % Copyright (C) 2017 Amir and Nili Beck
        %
        % This program is free software: you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation, either version 3 of the License, or
        % (at your option) any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program.  If not, see <http://www.gnu.org/licenses/>.

        saved_fmin_val = f(min_val) ;
        saved_fmax_val = f(max_val) ;

        if (saved_fmin_val * saved_fmax_val >0)
            error('f(lb)*f(ub)>0')
        end

        if (min_val > max_val)
            error ('minimal value is bigger than maximal_value\n') ;
        end

        iter=0;
        while (max_val-min_val>eps)

            %if the function is linear in this range, find the root
            r = (max_val-(saved_fmax_val/saved_fmin_val) * min_val) / (1-(saved_fmax_val/saved_fmin_val)) ;
            saved_froot = f(r) ;
            if (abs(saved_froot) < eps)
                out = r;
                return ;
            end

            iter=iter+1;
            mid=(min_val+max_val)/2;
            changed_limits = false ; %haven't changed the limits yet

            %checking if r maybe be one of the new range limits
            if (r > mid)
                %r may be the the new min_val
                if (saved_froot * saved_fmax_val <0 )
                    min_val = r ;
                    saved_fmin_val = saved_froot ;
                    changed_limits = true ;
                end
            else
                %r may be the new max_val
                if ( saved_froot * saved_fmin_val < 0)
                    max_val = r ;
                    saved_fmax_val = saved_froot ;
                    changed_limits = true ;
                end
            end

            if (~changed_limits)
                %didn't use r
                saved_fmid = f(mid) ;
                if(saved_fmin_val* saved_fmid >0)
                    min_val=mid;
                    saved_fmin_val = saved_fmid ;
                else
                    max_val=mid;
                    saved_fmax_val = saved_fmid ;
                end
            end
        end

        if (~changed_limits)
            out=mid ;
        else
            out = r;
        end

    end






end  % of function







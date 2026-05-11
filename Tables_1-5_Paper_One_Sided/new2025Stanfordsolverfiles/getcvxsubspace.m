function [d] = getcvxsubspace(X, Ds)
% Get a pre-conditioner that is the linear/affine combination of base pre-conditioners

X = X / trace(X);
k = size(Ds, 2); %#ok

cvx_begin sdp quiet
cvx_precision high

use_mosek=true;
try
    cvx_solver mosek
catch
    fprintf('mosek not available; using sdpt3 solver\n');
    cvx_solver sdpt3
    use_mosek=false;
end
variable b(k, 1)
variable tau nonnegative

if use_mosek
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 3600.0)
end

maximize tau
subject to
X - diag(Ds * b) >= 0; %#ok
diag(Ds * b) - X * tau >= 0; %#ok
cvx_end

d = Ds * b;

end % End function
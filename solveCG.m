function [u, niter, flag] = solveCG(funAx, f, s, tol, maxiter)
% SOLVECG   Conjugate Gradients method.
%
%    Input parameters: 
%           A : Symmetric, positive definite NxN matrix 
%           f : Right-hand side Nx1 column vector 
%           s : Nx1 start vector (the initial guess)
%         tol : relative residual error tolerance for break
%               condition 
%     maxiter : Maximum number of iterations to perform
%
%    Output parameters:
%           u : Nx1 solution vector
%       niter : Number of iterations performed
%        flag : 1 if convergence criteria specified by TOL could
%               not be fulfilled within the specified maximum
%               number of iterations, 0 otherwise (= iteration
%               successful).

% Author : Andreas Klimke, Universität Stuttgart
% Version: 1.0
% Date   : May 13, 2003
	
u = s;         % Set u_0 to the start vector s
As=feval(funAx,s);
r = f - As;   % Compute first residuum
p = r; 
r2=conj(r).*r;
rho = real(sum(sum(r2)));

niter = 0;     % Init counter for number of iterations
flag = 0;      % Init break flag

% Compute norm of right-hand side to take relative residuum as
% break condition.
% normf = norm(f);
% if normf < eps  % if the norm is very close to zero, take the
%                 % absolute residuum instead as break condition
%                 % ( norm(r) > tol ), since the relative
%                 % residuum will not work (division by zero).
%   warning(['norm(f) is very close to zero, taking absolute residuum' ... 
% 					 ' as break condition.']);
% 	normf = 1;
% end

while (sqrt(rho)> tol)   % Test break condition
	%a = A*p;
    a=feval(funAx,p);
    cap=conj(a).*p;
	alpha = rho/real(sum(sum(cap)));
	u = u + alpha*p;
	r = r - alpha*a;
    r2=conj(r).*r;
    rho_new = real(sum(sum(r2)));
	p = r + rho_new/rho * p;
	rho = rho_new;
	niter = niter + 1;
	if (niter == maxiter)         % if max. number of iterations
		flag = 1;                   % is reached, break.
		break
	end
end
function J = approximate_jacobian(fun, x)
% Approximate Jacobian of vector-valued function fun at point x

% INPUTS:
% fun: function handle, returns column vector
% x: column vector at which to evaluate Jacobian

% OUTPUTS:
% J: m x n Jacobian matrix

    delta_x = 1e-6; % small step
    f0 = fun(x); % evaluate function at original point
    m = length(f0); % number of outputs
    n = length(x); % number of variables
    J = zeros(m,n); % preallocate Jacobian
    
    e = zeros(n,1);

    for j = 1:n
        % Unit vector in j-th direction
        % e = zeros(n,1);
        e(j) = 1;

        % Evaluate function at x+h and x-h
        f_plus  = fun(x + delta_x*e);
        f_minus = fun(x - delta_x*e);

        % Central difference approximation for j-th column
        J(:,j) = (f_plus - f_minus) / (2*delta_x);

        e(j) = 0;
    end
end

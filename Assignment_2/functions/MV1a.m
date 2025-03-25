function [y, u] = MV1a_Controller(A, B, C, Ay, By, k, rho, omega, e, N)
    % MV1a_Controller: Implements MV1a control law iteratively.
    %
    % Inputs:
    %   A  - Coefficients of A(q^-1) polynomial (vector)
    %   B  - Coefficients of B(q^-1) polynomial (vector)
    %   C  - Coefficients of C(q^-1) polynomial (vector)
    %   Ay - Filter polynomial for y (vector)
    %   By - Filter polynomial for reference (vector)
    %   k  - Input delay (integer)
    %   rho - Regularization parameter (scalar)
    %   omega - Reference signal (vector of length N)
    %   e - Disturbance error signal (vector of length N)
    %   N  - Number of time steps (integer)
    %
    % Outputs:
    %   y - Output signal (vector of length N)
    %   u - Control input signal (vector of length N)

    % Solve the Diophantine equation for G and S
    [G, S] = diophantine(conv(Ay, A), conv(By, C), k);

    % Ensure input signals have length N
    omega = [omega, zeros(1, N - length(omega))]; 
    e = [e, zeros(1, N - length(e))];

    % Get polynomial orders
    na = length(A) - 1;
    nb = length(B) - 1;
    nc = length(C) - 1;
    ng = length(G) - 1;
    ns = length(S) - 1;

    % Initialize output and control input vectors
    y = zeros(1, N);
    u = zeros(1, N);

    % Extract first coefficient of B(q^-1) for regularization
    b0 = B(1);
    alpha = rho / b0;

    % Iterative MV1a control computation
    for t = max([na, nb, nc, ng, ns, k]) + 1 : N  
        % Compute y using past values (and include u in the equation)
        y_t = -sum(A(2:end) .* flip(y(t-na:t-1))) ...
              + sum(B(1:nb+1) .* flip(u(t-nb:t-1+k))) ...  % Adjusted here
              + sum(C .* flip(e(t-nc:t-1))) ...
            ;
        % Compute control input u using MV1a law
        u_t_1 = u(t-1); % H_u(q) = 1 - q^-1

        u(t) = (1 / (b0 + alpha)) * (... 
               -alpha * u_t_1 ...
               + sum(S .* flip(y(t-ns:t-1))) ...
               - sum(G .* flip(omega(t-ng:t-1))) ...
              );

        % Store output
        y(t-1) = y_t;
    end
end
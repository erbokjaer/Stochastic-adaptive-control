function [y, u] = MV1a(A, B, C, k, omega, e, N,  Ay, By, Au, Bu, Aw, Bw, rho, y_init, u_init)
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

    [G, S] = diophantine(A, 1, k);
    R1 = conv(Au,conv(B,G));
    R2 = rho*conv(C,Bu);

    % diff = numel(R1) - numel(R2);
    % if diff > 0
    %     R2 = [R2, zeros(1,abs(diff))];
    % elseif diff < 0 
    %     R1 = [R1, zeros(1,abs(diff))];
    % end

    R = R1;



    % Ensure input signals have length N
    omega = [y_init; omega]; 
    e = [zeros(numel(y_init),1) ; e];

    % Get polynomial orders
    na = numel(A) - 1;
    nb = numel(B) - 1;
    nc = numel(C) - 1;
    nr = numel(R) - 1;
    ng = numel(G) - 1;
    ns = numel(S) - 1;

    % Initialize output and control input vectors
    y = [y_init; zeros(N,1)];
    u = [u_init; zeros(N,1)];

    % Extract first coefficient of B(q^-1) for regularization
    b0 = B(1);
    alpha = rho / b0;

    % Iterative MV1a control computation
    for t = max([na, nb, nc, nr, ns, ng, ns, k]) + 1 : N + max(numel(y_init),numel(u_init))
        Y_t = y(t - (1:na));
        U_t = u(t - k - (0:nb));
        E_t = e(t - (0:nc));

        % Compute y using past values (and include u in the equation)

        y(t) = - A(2:end) * Y_t + B * U_t + C * E_t;

        Y_s = y(t - (0:ns));
        U_r = u(t - (1:nr));
        W_t = omega(t);

        % Compute control input u using MV1a law % H_u(q) = 1 - q^-1
        u(t) = 1/(R(1)) * (W_t-S*Y_s - R(2:end)*U_r);
        % Store output
    end

    y = y(numel(y_init)+1:end);
    u = u(numel(u_init)+1:end);
    
end
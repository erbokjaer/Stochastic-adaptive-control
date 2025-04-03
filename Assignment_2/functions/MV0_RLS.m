function [y,u] = MV0_RLS(A, B, omega, e, N, k, y_init, u_init, theta_init, P_init)
    C = 1;
    [G,S] = diophantine(A,C,k);
    R = conv(B,G);

    % Ensure input signals have length N
    omega = [y_init; omega]; 
    e = [zeros(numel(y_init),1) ; e];

    % Get polynomial orders
    na = numel(A) - 1;
    nb = numel(B) - 1;
    nr = numel(R) - 1;
    ng = numel(G) - 1;
    ns = numel(S) - 1;
    nc = numel(C) - 1;

    simStart = max([na, nb, nr, ns, ng, ns]) + k;

    % Initialize output and control input vectors
    y = [y_init; zeros(N,1)];
    u = [u_init; zeros(N,1)];

    num = na + nb + 1;
    % Initialize Recursive Least Squares (RLS) Parameters
    theta_rls = zeros(num, N); % Store estimated parameters over time
    P = P_init; % Initial covariance (high uncertainty)
    lambda = 1; % Forgetting factor

    theta_hat = theta_init; % Initial parameter estimates
    % theta_hat = [A(2:end),B]';
    P_plot = theta_hat;

    % Iterative MV1a control computation
    for t = simStart: N + max(numel(y_init),numel(u_init))
        Y_t = y(t - (1:na));
        U_t = u(t - k - (0:nb));
        E_t = e(t - (0:nc));

        % Compute y using past values (and include u in the equation)
        y(t) = - A(2:end) * Y_t + B * U_t + C * E_t;

        Phi_t = [-Y_t', U_t']';

        % Kalman gain
        K_t = P * Phi_t / (lambda + Phi_t' * P * Phi_t);
        
        % Prediction error
        e_t = y(t) - Phi_t' * theta_hat;
        
        % Parameter update
        theta_hat = theta_hat + K_t * e_t;
        
        % Covariance update
        P = (P - K_t * Phi_t' * P) / lambda;
        P_plot(:, t) = diag(P);
        
        
        % Store estimates over time
        theta_rls(:, t) = theta_hat;
        A_est = [1,theta_hat(1:na)'];
        B_est = theta_hat(na + 1:end)';

        [G,S] = diophantine(A_est,C,k);
        R = conv(B_est,G);
        

        Y_s = y(t - (0:ns));
        U_r = u(t - (1:nr));
        W_t = omega(t);
        % W_t = omega(t);
        
        % Compute control input u using MV1a law % H_u(q) = 1 - q^-1
        u(t) = 1/(R(1)) * (W_t - S * Y_s - R(2:end) * U_r);
        % u(t) = 0;

        % Store output
    end
    y = y(numel(y_init)+1:end);
    u = u(numel(u_init)+1:end);

    disp("A_est_final = [" + num2str(A_est) + "]")
    disp("B_est_final = [" + num2str(B_est) + "]")
end
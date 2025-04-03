function [y, theta_rls] = armax_forgetting(A, B, C, u, e, N, y_init, theta_init, P_init, k, forgetting_method, lambda)
   

    % Get polynomial orders
    na = numel(A) - 1;
    nb = numel(B) - 1;
    nc = numel(C) - 1;

    simStart = max([na, nb + k, nc]) + 1;

    % Initialize output vector
    y = [y_init; zeros(N, 1)];

     % Ensure input signals have length N
    u = [zeros(numel(y_init) + k, 1); u]; 
    e = [zeros(numel(y_init), 1); e];

    % Initialize Recursive Least Squares (RLS) parameters
    num_params = na + (nb + 1) + nc;
    theta_rls = zeros(num_params, N); % Store estimated parameters over time
    theta_hat = theta_init; % Initial parameter estimates
    P = P_init; % Initial covariance (high uncertainty)

    % Iterative parameter estimation
    for t = simStart:N + numel(y_init)
        % Construct regressor vector
        Y_t = y(t - (1:na));
        U_t = u(t - k - (0:nb));
        E_t = e(t - (1:nc));

        Phi_t = [-Y_t', U_t', E_t']';  % Regressor vector

        % Output computation from current parameters (simulation model)
        y(t) = - A(2:end) * Y_t + B * U_t + C * [e(t);E_t];

        % Prediction error
        e_t = y(t) - Phi_t' * theta_hat;
        
        
        % Covariance update based on forgetting method
        if strcmp(forgetting_method, 'exponential')
            % Kalman gain
            K_t = P * Phi_t / (lambda + Phi_t' * P * Phi_t);
            
            % Parameter update
            theta_hat = theta_hat + K_t * e_t;
            
            P = (P - K_t * Phi_t' * P) / lambda;
        
        else
            % Kalman gain
            K_t = P * Phi_t / (1 + Phi_t' * P * Phi_t);
            
            % Parameter update
            theta_hat = theta_hat + K_t * e_t;
            
            P = P - K_t * Phi_t' * P; % Standard RLS update
        end
        
        % Store estimates over time
        theta_rls(:, t) = theta_hat;


        
        
        

    end

    % Remove initial conditions from output
    y = y(numel(y_init) + 1:end);
    theta_rls = theta_rls(:, numel(y_init) + 1:end);
end

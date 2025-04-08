function [y, u, theta_rls] = MVa1_RLS(A, B, C, k, omega, e, N,theta_init, P_init,  Ay, By, Au, Bu, Aw, Bw, rho, y_init, u_init, forgetting_method, lambda)
    

    % Solve the Diophantine equation for G and S
    [G, S] = diophantine(conv(Ay, A), conv(By, C), k);

    % [G, S] = diophantine(A, 1, k);
    b0 = B(1);
    alpha = rho / b0;
    R1 = conv(Au,conv(B,G));
    R2 = alpha*conv(C,Bu);
    
    
    diff = numel(R1) - numel(R2);
    if diff > 0
        R2 = [R2, zeros(1,abs(diff))];
    elseif diff < 0 
        R1 = [R1, zeros(1,abs(diff))];
    end

    R = conv(Aw,conv(Ay,R1 + R2));
    % R1_conv = conv(conv(Ay,R1),Aw);
    % R2_conv = conv(conv(Ay,R2),Aw);
    % R = R1_conv + R2_conv;

    W_in = conv(conv(conv(Au,Ay),Bw),C);
    Y_in = conv(Au,conv(Aw,S));


    % Ensure input signals have length N
    omega = [y_init; omega]; 
    e = [zeros(numel(y_init),1) ; e];
    E_t_est = zeros(numel(e),1);

    % Get polynomial orders
    na = numel(A) - 1;
    nb = numel(B) - 1;
    nc = numel(C) - 1;
    nr = numel(R) - 1;
    ng = numel(G) - 1;
    ns = numel(S) - 1;
    nw_in = numel(W_in) - 1;
    ny_in = numel(Y_in) - 1;
    simStart = max([na, nb, nc, nr, ns, ng, ns,nw_in,ny_in]) + k;

    % Initialize output and control input vectors
    y = [y_init; zeros(N,1)];
    u = [u_init; zeros(N,1)];

     % Initialize Recursive Least Squares (RLS) parameters
    num_params = na + (nb + 1) + nc;
    theta_rls = zeros(num_params, N); % Store estimated parameters over time
    theta_hat = theta_init; % Initial parameter estimates
    P = P_init; % Initial covariance (high uncertainty)


    % Iterative MV1a control computation
    for t = simStart: N + max(numel(y_init),numel(u_init))
        Y_t = y(t - (1:na));
        U_t = u(t - k - (0:nb));
        E_t = e(t - (0:nc));

        % Compute y using past values (and include u in the equation)
        y(t) = - A(2:end) * Y_t + B * U_t + C * E_t;
        
        % Construct regressor vector
        Phi_t = [-Y_t', U_t', E_t_est(t - (1:nc))']';  % Regressor vector
        
        % Prediction error
        e_t = y(t) - Phi_t' * theta_hat;

        E_t_est(t) = e_t;
        
        
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

        % if mod(t,50)
        %     P = P_init;
        % end 
        
        % Store estimates over time
        theta_rls(:, t) = theta_hat;
        
        A_est = [1, theta_hat(1:na)'];
        B_est = theta_hat(na+1:na+1+nb)';
        C_est = [1, theta_hat(na+1+nb+1:end)'];


            % Solve the Diophantine equation for G and S
        [G, S] = diophantine(conv(Ay, A_est), conv(By, C_est), k);
    
        % [G, S] = diophantine(A, 1, k);
        b0 = B_est(1);
        alpha = rho / b0;
        R1 = conv(Au,conv(B_est,G));
        R2 = alpha*conv(C_est,Bu);
        
        
        diff = numel(R1) - numel(R2);
        if diff > 0
            R2 = [R2, zeros(1,abs(diff))];
        elseif diff < 0 
            R1 = [R1, zeros(1,abs(diff))];
        end
    
        R = conv(Aw,conv(Ay,R1 + R2));
        W_in = conv(conv(conv(Au,Ay),Bw),C_est);
        Y_in = conv(Au,conv(Aw,S));
        
        Y_s = y(t - (0:ny_in));
        U_r = u(t - (1:nr));
        W_t = omega(t - (0:nw_in));
        % W_t = omega(t);
        
        % Compute control input u using MV1a law % H_u(q) = 1 - q^-1
        u(t) = 1/(R(1)) * (W_in * W_t - Y_in * Y_s - R(2:end) * U_r);

        % Store output
    end

    y = y(numel(y_init)+1:end);
    u = u(numel(u_init)+1:end);
end
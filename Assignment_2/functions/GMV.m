function [y, u] = GMV(A, B, C, k, omega, e, N,  Ay, By, Au, Bu, Aw, Bw, rho, y_init, u_init)
    

    % Solve the Diophantine equation for G and S
    [G, S_diop] = diophantine(conv(Ay, A), conv(By, C), k);
    disp("G:")
    disp(G)

    disp("S:")
    disp(S_diop)


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

    Q = conv(conv(conv(Au,Ay),Bw),C);
    S = conv(Au,conv(Aw,S_diop));


    % Ensure input signals have length N
    omega = [y_init; omega]; 
    e = [zeros(numel(y_init),1) ; e];

    % Get polynomial orders
    na = numel(A) - 1;
    nb = numel(B) - 1;
    nc = numel(C) - 1;
    nr = numel(R) - 1;
    ng = numel(G) - 1;
    ns = numel(S_diop) - 1;
    nw_in = numel(Q) - 1;
    nS_in = numel(S) - 1;
    simStart = max([na, nb, nc, nr, ns, ng, ns,nw_in,nS_in]) + k;

    % Initialize output and control input vectors
    y = [y_init; zeros(N,1)];
    u = [u_init; zeros(N,1)];

    % Iterative MV1a control computation
    for t = simStart: N + max(numel(y_init),numel(u_init))
        Y_t = y(t - (1:na));
        U_t = u(t - k - (0:nb));
        E_t = e(t - (0:nc));

        % Compute y using past values (and include u in the equation)
        y(t) = - A(2:end) * Y_t + B * U_t + C * E_t;

        Y_s = y(t - (0:nS_in));
        U_r = u(t - (1:nr));
        W_t = omega(t - (0:nw_in));
        % W_t = omega(t);
        
        % Compute control input u using MV1a law % H_u(q) = 1 - q^-1
        u(t) = 1/(R(1)) * (Q * W_t - S * Y_s - R(2:end) * U_r);

        % Store output
    end
    y = y(numel(y_init)+1:end);
    u = u(numel(u_init)+1:end);
end
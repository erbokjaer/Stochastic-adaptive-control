function struct = run_case(s, SP, Stochastic, ctrl)
    struct = s;
    A = s.A;
    B = s.B;
    N = s.N;
    k = s.k;
    y_init = s.y_init;
    u_init = s.u_init;

    switch SP
        case 1
            omega = zeros(N,1); % Constant reference setpoint
        case 2
            omega = ones(N,1);
        case 3
            T = N/2; % Periode in sampels
            omega = square(2*pi*1/T*(1:N),50)';
        case 4
            omega = s.ref;
        otherwise
            msg = "Choose a valid value for: Setpoint";
            error(msg)
    end
    struct.omega = omega;
    
    switch Stochastic
        case 1
            e = s.sigma_e * randn(N, 1);  % Process noise (Gaussian)
        case 2
            e = zeros(N,1);
        otherwise
            msg = "Choose a valid value for: Stochastic";
            error(msg)
    end
    struct.e = e;
    
    switch ctrl
        case 1
            [y,u] = Pctrl(A, B, s.alpha, omega, e, N, k);
        case 2
            [y,u] = MV(A, B, e, N, k);
        case 3
            [y,u] = MV0(A, B, omega, e, N, k, y_init, u_init);
        case 4
            [y,u] = ARX(A, B, e, N, k, y_init, u_init, s.u_data);
        case 5
            [y,u] = ARX(A, B, e, N, k, y_init, u_init, zeros(N,1));
        otherwise
            msg = "Choose a valid value for: Control type";
            error(msg)
    end
    struct.y = y;
    struct.u = u;
end


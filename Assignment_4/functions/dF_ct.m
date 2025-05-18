function JF = dF_ct(x0)
    % Jacobian of F_ct w.r.t x
    kappa = 5*10^-5;
    lambda = 3;
    Lambda = 5*10^-5;
    beta = 0.0065;
    H = 0.05;
    JF = [
        (x0(3) - beta)/Lambda,  lambda, x0(1)/Lambda; % Corrected (3,3) entry
        beta/Lambda, - lambda, 0;
        -kappa*H, 0, 0;
    ];
end
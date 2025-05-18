function dxdt = F_ct(x0)
    kappa = 5*10^-5;
    lambda = 3;
    Lambda = 5*10^-5;
    beta = 0.0065;
    H = 0.05;
    dxdt = [
        (x0(3) - beta)/Lambda * x0(1) + lambda * x0(2);
        beta/Lambda * x0(1) - lambda * x0(2);
        -kappa*H*x0(1);
    ];
end




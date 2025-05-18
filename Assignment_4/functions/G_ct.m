function B_ct = G_ct(x0)
    Lambda = 5*10^-5;
    B_ct = [
        x0(1)/Lambda; % Corrected (1,1) entry from G
        0;
        0;
    ];
end


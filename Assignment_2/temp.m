omega_filt = 1
G = tf([omega_filt],[0.1,omega_filt]);
H = c2d(G,0.1);
input = H.Numerator{1}; 
% Aw = H.Denominator{1};
% Bw = input(find(input~=0,1,'first'):find(input~=0,1,'last'));
% Ay = H.Denominator{1};
% By = input(find(input~=0,1,'first'):find(input~=0,1,'last'));





% [s24.y,s24.u] = ARMAX(s24.A, s24.B, s24.C, s24.k,s24.e, s24.N, s24.y_init, s24.u_init);

% [s11.y,s11.u] = ARX(s11.A, s11.B, s11.e, s11.N, s11.k, s11.y_init, s11.u_init, zeros(s11.N,1));
% [s11.y,s11.u] = MV1a(s11.A, s11.B, s11.C, 1, 1, s11.k, s11.rho , s11.omega,s11.e, s11.N);

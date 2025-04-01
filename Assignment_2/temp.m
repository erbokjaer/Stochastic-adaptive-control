omega_filt = 1
G = tf([omega_filt],[0.1,omega_filt]);
H = c2d(G,0.1);
input = H.Numerator{1}; 
% Aw = H.Denominator{1};
% Bw = input(find(input~=0,1,'first'):find(input~=0,1,'last'));
% Ay = H.Denominator{1};
% By = input(find(input~=0,1,'first'):find(input~=0,1,'last'));

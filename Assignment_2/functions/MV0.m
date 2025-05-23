function [y,u] = MV0(A, B, omega, e, N, k, y_init, u_init)
        [G,S] = diophantine(A,1,k);

        R = conv(B,G);
        e = [zeros(max([numel(A),numel(B)]),1);e];
        u = [u_init; zeros(N+max([numel(A),numel(B)]),1)];
        y = [y_init; zeros(N+max([numel(A),numel(B)]),1)];
        omega = [zeros(max([numel(A),numel(B)]),1);omega];

    for t = max([length(A), length(B),length(R)])+k:N + max([numel(A),numel(B),numel(R)])
        % Collect past values safely
        y_past = y(t-1:-1:t-length(A)+1);
        u_past = u(t-k:-1:t-k-length(B)+1);

        % Compute system output
        y(t) = -A(2:end) * y_past + B * u_past + e(t);


        % MVC0
        u(t) = 1/(R(1)) * (omega(t)-S*y(t:-1:t-length(S)+1) - R(2:end)*u(t-1:-1:t-length(R)+1));
    end
    y = y(numel(y_init):end - 1 - max([numel(A),numel(B),numel(R)]));
    u = u(numel(u_init):end - 1 - max([numel(A),numel(B),numel(R)]));

    % Alternative implimentation throug GMV
    % Ay = [1];
    % By = [1]; 
    % Au = [1]; 
    % Bu = [1]; 
    % Aw = [1]; 
    % Bw = [1]; 
    % C = 1;
    % rho = 0; 
    % 
    % 
    % [y, u] = GMV(A, B, C, k, omega, e, N, Ay, By, Au, Bu, Aw, Bw, rho, y_init, u_init);

    
end


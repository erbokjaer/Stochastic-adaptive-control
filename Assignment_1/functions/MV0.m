function [y,u] = MV0(A, B, omega, e, N, k)
        [G,S] = diophantine(A,1,k);
        
        R = conv(B,G);

        u = zeros(N,1);
        y = zeros(N,1);
    for t = max([length(A), length(B),length(R)])+1:N
        % Collect past values safely
        y_past = y(t-1:-1:t-length(A)+1);
        u_past = u(t-k:-1:t-k-length(B)+1);
    
        % Compute system output
        y(t) = -A(2:end) * y_past + B * u_past + e(t);
        
        % MVC0
        u(t) = 1/(R(1)) * (omega(t)-S*y(t:-1:t-length(S)+1) - R(2:end)*u(t-1:-1:t-length(R)+1));
    end
end


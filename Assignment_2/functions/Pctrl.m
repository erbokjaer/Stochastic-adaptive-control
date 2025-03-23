function [y,u] = Pctrl(A, B, alpha, omega, e, N, k)
    u = zeros(N,1);
    y = zeros(N,1);

    for t = max([length(A), length(B)])+1:N
        % Collect past values safely
        y_past = y(t-1:-1:t-length(A)+1);
        u_past = u(t-k:-1:t-k-length(B)+1);
    
        % Compute system output
        y(t) = -A(2:end) * y_past + B * u_past + e(t);
        
        % Compute control input using the P-controller
        u(t) = alpha * (omega(t) - y(t));
    end
end

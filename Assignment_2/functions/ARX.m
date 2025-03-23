function [y, u] = ARX(A, B, e, N, k, y_init, u_init, u_in)
    u = [u_init; u_in(1:N-numel(u_init))];
    y = [y_init; zeros(N-numel(y_init),1)];

    for t = max([length(A), length(B)])+k:N
        % Collect past values safely
        y_past = y(t-1:-1:t-length(A)+1);
        u_past = u(t-k:-1:t-k-length(B)+1);
    
        % Compute system output
        y(t) = -A(2:end) * y_past + B * u_past + e(t);
    end
end



function [y, u] = ARMAX(A, B, C, k, e, N, y_init, u_init)

    e = [zeros(numel(y_init),1) ; e];

    % Get polynomial orders
    na = numel(A) - 1;
    nb = numel(B) - 1;
    nc = numel(C) - 1;
    simStart = max([na, nb, nc]) + k;

    % Initialize output and control input vectors
    y = [y_init; zeros(N,1)];
    u = [u_init; zeros(N,1)];
    e = [zeros(numel(y_init),1) ; e];

    % Iterative MV1a control computation
    for t = simStart: N + max(numel(y_init),numel(u_init))
        Y_t = y(t - (1:na));
        U_t = u(t - k - (0:nb));
        E_t = e(t - (0:nc));

        % Compute y using past values (and include u in the equation)

        y(t) = - A(2:end) * Y_t + B * U_t + C * E_t;

    end

    y = y(numel(y_init)+1:end);
    u = u(numel(u_init)+1:end);
end


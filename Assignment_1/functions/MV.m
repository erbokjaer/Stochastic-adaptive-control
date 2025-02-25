function [y,u] = MV(A, B, e, N, k)
    [y,u] = MV0(A, B, zeros(N,1), e, N, k);
end


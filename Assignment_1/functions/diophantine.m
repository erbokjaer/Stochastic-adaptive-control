function [G, S] = diophantine(A,B,m)
    G = [];
    S = [B, zeros(1,length(A)-length(B))]; % Pad B with zeros to make S as long as A
    for i = 1:m 
        % Augment with first element of S
        G = [G, S(1)];
    
        % Update S
        S = [S(2:end) - S(1) *A(2:end), 0] ;
    end
    
    % Remove last element
    S = S ( 1 : end-1);
end
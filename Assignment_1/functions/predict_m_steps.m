function y_pred = predict_m_steps(A, B, y, u_data, m, N, k)
    [G,S] = diophantine(A,1,m);
    R = conv(B,G);

    % Initialize predicted values as zeros before sim
    offset = numel(B)+max([numel(S),numel(B),numel(R)])+k;
    
    % Store predictions
    y_pred = [zeros(1, offset),zeros(1, N)];  
    
    % Use simulation data (s.y) for prediction:
    y_sampl =  [zeros(offset,1);y];
    
    % Use same 
    u_sampl = [zeros(offset,1);u_data];
    
    % Simulate predictions iteratively
    for t = numel(B)+max([numel(S),numel(B),numel(R),offset + numel(B)])+k:N + 2 * numel(B) + offset
        y_pred(t+m-offset) =  R * u_sampl(t+m-numel(B)-k-offset-(0:numel(R)-1)) + S * y_sampl(t-offset-(0:numel(S)-1));
    end
    y_pred = y_pred(offset+1:N+ offset)';
end


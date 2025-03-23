function y_pred = predict_m_steps(A, B, y, u_data, m, N, k)
    [G,S] = diophantine(A,1,m);
    R = conv(B,G);
    nb = numel(B)-1;

    % Initialize predicted values as zeros before sim
    offset = numel(B)+max([numel(S),numel(B),numel(R)])+k;
    
    % Store predictions
    y_pred = [zeros(1, offset),zeros(1, N)];  
    
    % Use simulation data (s.y) for prediction:
    y_sampl =  [zeros(offset,1);y];
    
    % Use same 
    u_sampl = [zeros(offset,1);u_data];

    end_sim = N+ 2*numel(B) + offset;
    start_sim = 2*numel(B)+offset+k;
    
    % Simulate predictions iteratively
    for t = start_sim:end_sim
        R_term = R * u_sampl(t+m-nb-k-offset-(1:numel(R)));
        S_term = S * y_sampl(t-offset-(0:numel(S)-1));
        y_pred(t+m-offset) =  R_term + S_term;
    end


    y_pred = y_pred(offset + 1:N + offset)';
end
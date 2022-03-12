function [c_out, c_outeq] = c_constrain(z)
    ALPHA = 0.2;
    BETA = 20;
    LAMBDA_T = (2*pi) / 3;
    N_SAMPLES = 40;
    N_STATES = 6;
    
    % z = (x1 ..... xn, u0, .... un-1)
    
    % travel_pos = 1;
    % elevation_pos = 5;
    
    % Q = 320x320
    % Q1s 1:240
    % R1s 241:320 
    
    travels = z(1 : N_STATES : N_SAMPLES*N_STATES); % 1 : 6 : 40*6 = 1:6:240
    elevatn = z(5 : N_STATES : N_SAMPLES*N_STATES); % 5 : 6 : 40*6 = 5:6:240
    
    % construct c
    c_out = zeros(N_SAMPLES, 1);
    c_outeq = [];
    
    for i = 1 : N_SAMPLES
        c_out(i) = ALPHA * exp( -BETA*(travels(i) - LAMBDA_T)^2) - elevatn(i);
    end
end


function options = create_options(varargin)
    % 
    % Returns
    % options.pi_0
    % options.tau_0
    % options.est_pi0
    
    options.pi_0 = 0;  % Assume half of all findings are null by default.
    options.tau_0 = .2; % T statistic threshold for null hypothesis
    options.est_pi0 = [];
    
    
end
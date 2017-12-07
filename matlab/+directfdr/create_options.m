function options = create_options(varargin)
    % 
    % Returns
    % options.pi_0
    % options.tau_0
    % options.est_pi0
    % options.nBoot
    
    options.pi_0 = 0;  % Rejection region for p-values
    options.tau_0 = .2; % Conservative non-rejection region for test statistics
    options.est_pi0 = [];
    
    % Display on or off; 
    options.verbose = 0;
    
    % Bootstrap options for parameter tuning
    options.nBoot = 100;
    options.Tau = linspace(.05,.5,4);
    
    % Speed up. 
    options.topk = 0; % Set to 5 or 10 to speed up bootstrap iterations;
    
end
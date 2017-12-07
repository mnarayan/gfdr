function [results, opts] = optimize(Tobs,Tperm,varargin)
    
    
    opts = directfdr.create_options();
    opts.verbose = 1;
    fdr_all = zeros(length(opts.Tau));
    [fdr_all results_all] = compute_fdr_tau(Tobs,Tperm,opts);
    min_fdr = min(fdr_all);
    
    % Perform bootstrapped pFDR estimation.
    m = length(Tobs);
    currstate = rng;
    fdr_boot = zeros(opts.nBoot,length(opts.Tau));
    opts.topk = 5; % For speed up. 
    opts.verbose = 0;
    for bootNo=1:opts.nBoot
        boot_idx = randsample(m,m,1); 
        fdr_boot(bootNo,:) = compute_fdr_tau(...
                                        Tobs(boot_idx), ...
                                        Tperm(boot_idx,:), ...
                                        opts);
    end
    
    
    % MSE Scores and Min. MSE results
    mse_score = compute_mse(fdr_boot,min_fdr);   
    [best_mse, best_idx] = min(mse_score);
    tau_best = opts.Tau(best_idx);
    
    disp('Optimized Complement of Rejection Region')
    disp(sprintf('(0, %.3f)',tau_best));
    
    results = results_all{best_idx};
    
    opts.bootstate = currstate;
    opts.mse_score = mse_score;
    opts.fdr_boot = fdr_boot; 
    opts.min_fdr = min_fdr;
    opts.tau_best = tau_best;
    opts.fdr_best = results;
    
end


function mse_scores = compute_mse(fdr_boot,min_fdr)
    % 
    % Inputs
    % fdr_boot is nBoot x nTau matrix of fdr values
    % min_fdr is the min. fdr value across all nTau 
    
    mse_scores = mean((fdr_boot - min_fdr).^2,1); 
    
end


function min_fdr = compute_min_fdr(fdr)
    
    min_fdr = min(fdr);
     
end


function [fdr_all varargout] = compute_fdr_tau(Tobs,Tperm,opts)
   
    fdr_type = 'est_fdr';
    for tauNo=1:length(opts.Tau)
        opts.tau_0 = opts.Tau(tauNo); 
        if(nargout>=2)
            [~,~,results{tauNo}] = directfdr.run(Tobs,Tperm,opts);
            tmp_fdr = results{tauNo}.(fdr_type);
            nonzero = find(tmp_fdr~=0);
            fdr_all(tauNo) = min(tmp_fdr(nonzero));
        else
            [~,~,results] = directfdr.run(Tobs,Tperm,opts);
            tmp_fdr = results.(fdr_type);
            nonzero = find(tmp_fdr~=0);
            fdr_all(tauNo) = min(tmp_fdr(nonzero));
        end
    end  
    
    if(nargout>=2)
        varargout{1} = results;
    end
    
end

function compute_ci(T)
    

    
end
function [fdr_hat pval_adj results opts] = run(t,tk,varargin)
    % Estimate False Discovery Rate using Storey's permutation method
    % 
    % 
    % Input
    %   t  - observed test statistic of length n_hyp x 1
    %   tk - permuted test statistics of length n_hyp x n_perm
    %   opts - (optional) defaults to create_options()
    % 
    % Output
    %   fdr_hat     - Estimated FDR
    %   pval_adj    - Adjusted p-values
    %   results     - Tabulated results
    %   opts        - options structure with intermediate parameteres;
    % %
    %
    %     
    % Description: 
    % 
    % Notation Table: 
    %               | Accept | Reject | Total
    %   Null True   |   U    |   V    | m0
    %   Alt True    |   T    |   S    | m1
    %               |   W    |   R    | m 
    %
    % Conservative estimate of proportion of null hypotheses: 
    % est_pi0(lambda) =             W(lambda)
    %                    ----------------------------------
    %                        (1-lambda) * (# of hypotheses) 
    % where W(lambda) = #{p_i > lambda}, 
    % 
    % Here lambda is threshold that determines the observe rejection region of interest, so that null p-values and belongs to (lambda,1). 
    % Alternatively, the above can also be indexed in terms of tau_0 threshold to formulate rejection regions via test statistics 
    % 
    % 
    %                est_pi0(lambda) * gamma        W(lambda) * gamma
    % Est_FDR       = -----------------------  =  ----------------------------
    %                     Pr(P<=gamma)          (1-lambda) * (R(gamma) V 1 )
    % 
    % where R(gamma) = #{p_i <= gamma}. Analogously, we use the parameter c to denote the threshold that determines rejection region for test statistics instead of the p-values
    % 
    % The implementation in this function follows Storey & Tibshirani (2001) using general test statistics rather than the p-value method of Storey (2002). 
    % 
    % 
    % Reference
    % Storey, J. D., & Tibshirani, R. (2001). Estimating the positive false discovery rate under dependence, with applications to DNA microarrays. Technical Report 2001-28, Department of Statistics, Stanford University, 2001.
    % 
    % Storey, J. D. (2002), A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64: 479–498. doi:10.1111/1467
    % 
    % Copyright 2015, 2017 Manjari Narayan
    % BSD-3 License
    
    if(nargin==2)
        opts = directfdr.create_options();
    elseif(nargin>=3)
        opts = varargin{1};
    end
    
    c = opts.pi_0;
    tau_0 = opts.tau_0;
    
	[tau,sorted_idx] = sort(abs(t),'descend');
	if(c==0)
		for ii=1:length(tau)
			R_obs(ii) = max(1,length(find(abs(t)>=tau(ii))));
		end
	else
		R_obs = max(1,length(find(abs(t)>=c)));
	end		
	W_obs = length(find(abs(t)<=tau_0));	
	
	if(c==0)
		R_hat = zeros(length(t),1);
		fdr_hat = zeros(length(t),1);	
		pval_adj = fdr_hat;
		W_hat = 0;
		val_0 = (abs(tk)<tau_0);
        % sum(val_0(:,1:10)==1,1)
        % size(sum(val_0(:,1:10,1)))
		W_hat = mean(sum(val_0==1,1));
		est_pi0 = W_obs./W_hat;
        if(opts.verbose)
            disp(sprintf('Est Pi_0: %.2f',est_pi0));
        end
        if(opts.topk)
            n_tau = opts.topk;
        else
            n_tau = length(tau); 
        end
        
		for ii=1:n_tau
			[val] = (abs(tk)>tau(ii));
			R_hat(ii) = mean(sum(val==1,1));
            fdr_hat(ii) = (R_hat(ii)./R_obs(ii)) * est_pi0;
            if(opts.verbose)            
                if(ii==1)
                %     disp('R_hat'); R_hat(ii)
                %     disp('R_obs'); R_obs(ii)
                %     W_obs
                %     W_hat
                     disp(sprintf('Smallest Est. FDR: %.7f',fdr_hat(ii)));
                end
            end
		end	
	
		pval_adj(end) = fdr_hat(end);
		for ii=1:(length(fdr_hat)-1)
			pval_adj(end-ii) = min(fdr_hat(end-ii),pval_adj(end-ii+1));
		end
	else
		pval_adj = NaN;
        est_pi0 = NaN;
		V_hat = sum(sum(abs(tk)>c))/size(tk,2); % does not use estimate of pi_0
		fdr_hat = V_hat/R_obs;
	end

    results = table();
    results.id = sorted_idx;
    results.est_fdr = fdr_hat;
    results.pval_adj = pval_adj;
    results = sortrows(results,{'id'});
    
    fdr_hat(sorted_idx) = fdr_hat;
    pval_adj(sorted_idx) = pval_adj;
    
    opts.est_pi0 = est_pi0;

end
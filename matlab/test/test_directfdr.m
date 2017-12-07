function test_directfdr
    
    
    % Simulate data
    m = 1000;
    m0 = .50*m;  
    n = 50;
    n1 = n/2+1;
    mu = .5;
    
    rng('default');
    X = randn(m,n); 
    X(m0:end,n1:end) = X(m0:end,n1:end) + mu;
    
    [~,~,~,stats] = ttest2(X(:,1:n1-1)',X(:,n1:end)');
    Tobs = stats.tstat;
    Tobs = reshape(Tobs,[m 1]);
    
    n_perm = 2000; 
    Tperm = zeros(m,n_perm); 
    for perm=1:n_perm
        perm_idx = randperm(n,n); 
        [~,~,~,stats] = ttest2(X(:,perm_idx(1:n1-1))',X(:,perm_idx(n1:end))');
        Tperm(:,perm) = stats.tstat;
    end
    
    opts = directfdr.create_options();
    [est_fdr pval_adj results] = directfdr.run(Tobs,Tperm,opts);
    
    disp('No of rejected tests ...')
    sum(pval_adj<=.05)
    
    disp('Results with Est_FDR<.01')    
    results(find(results.est_fdr<.01),:)
    
    % est_fdr(results.est_fdr<.01)
    % pval_adj(results.est_fdr<.01)
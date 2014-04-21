% This script tests the PSVD solvers on the 'tdmat_news20' matrix with 
% a varying number of singular triplets

addpath('/users/zwang/research/psvd_test/LMSVD-root')     % LMSVD
addpath('/users/zwang/research/psvd_test/propack_PSVD')   % PROPACK
addpath('/users/zwang/research/psvd_test/sparseSVD-root') % sparse SVD
addpath('/users/zwang/research/psvd_test/bchdav_psvd')    % bchdav
warning off

% load the matrix for PSVD
load('~/research/data/tdmat_news20.mat') % tdmat_news20: 53,975 x 11,269

% increase the number of singular triplets to compute and record
% the CPUtime for each method
kvals = [400 800 1200 1600 2000];
% compute the F-norm of the matrix 
data_fnorm = norm(data, 'fro');
% compute the estimated 2-norm of the matrix
data_2norm = normest(data);
% specify the convergence tolerance if required by the algorithm
tol = 1e-8; 
opts = struct('tol',tol);

%======================================================================
% randomized SVD (power method)
% oversampling param = 2, power param = 2
fprintf('rsvd_p with p=2 q=2\n')
fprintf('k    time    err_matrix    err_vector_max\n')
for i = 1:length(kvals)
    k = kvals(i); % k: number of singular triplets
    %--------------------------------------------------
    ts = cputime();
    [U,S,V] = rsvd_p(data, k, 2, 2, 'std_normal');
    te = cputime()-ts;
    %--------------------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*S*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*S;
    err_vec = zeros(k,1);
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    % output for analysis
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end

%======================================================================
% randomized SVD (power method)
% oversampling param = k, power param = 2
fprintf('rsvd_p with p=k q=2\n')
fprintf('k    time    err_matrix    err_vector_max\n')
for i = 1:length(kvals)
    k = kvals(i); % k: number of singular triplets
    err_vec = zeros(k,1);
    %---------------------------------------------------
    ts = cputime();
    [U,S,V] = rsvd_p(data, k, k, 2, 'std_normal');
    te = cputime()-ts;
    %---------------------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*S*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*S;
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    % output for analysis
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end

%======================================================================
% randomized SVD (Krylov subspace)
% oversampling param = 2, power param = 2 (default setting in pca.m)
fprintf('rsvd_k with p=2 q=2\n')
fprintf('k    time    err_matrix    err_vector_max\n')
for i = 1:length(kvals)
    k = kvals(i); % k: number of singular triplets
    err_vec = zeros(k,1);
    %--------------------------------------------------
    ts = cputime();
    [U,S,V] = rsvd_k(data, k, 2, 2, 'std_normal');
    te = cputime()-ts;
    %--------------------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*S*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*S;
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end

%======================================================================
% LBD with partial reorthogonalization (PROPACK)
fprintf('lansvd (PROPACK): no implicit restart \n')
fprintf('k    time    err_matrix    err_vector_max\n')
for i = 1:length(kvals)
    k = kvals(i); % k: number of singular triplets
    err_vec = zeros(k,1);
    %----------------------------------------
    ts = cputime();
    [U,S,V] = lansvd(data, k, 'L', opts);
    te = cputime()-ts;
    %----------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*S*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*S;
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end

%======================================================================
% built-in function svds in MATLAB (eigs on [0 A; A' 0])
fprintf('svds (MATLAB built-in function) \n')
fprintf('k    time    err_matrix    err_vector_max\n')
for i = 1:length(kvals)
    k = kvals(i); % k: number of singular triplets
    err_vec = zeros(k,1);
    %----------------------------------------
    ts = cputime();
    [U,S,V] = svds(data,k,'L',opts);
    te = cputime()-ts;
    %----------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*S*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*S;
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end

%======================================================================
% limited memory svd 
fprintf('limited memory svd \n')
fprintf('k    time    err_matrix    err_vector_max\n')
for i = 1:length(kvals)
    k = kvals(i); % k: number of singular triplets
    err_vec = zeros(k,1);
    %-------------------------------------
    ts = cputime();
    [U,S,V] = lmsvds(data,k,opts);
    te = cputime()-ts;
    %-------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*S*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*S;
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end

%======================================================================
% bchdav 
fprintf('bchdav \n')
fprintf('k    time    err_matrix    err_vector_max\n')
opts_b = struct('polym',6,'blk',50,'vimax',250);
for i = 1:length(kvals)
    
    k = kvals(i); % k: number of singular triplets
    %--------------------------------------------------
    ts = cputime();
    [U, S, V] = bchdav_psvd(data, k, opts_b);
    te = cputime()-ts;
    %--------------------------------------------------
    % compute the low-rank approximation error in F-norm
    err_mat = norm(data - U*diag(S)*V','fro')/data_fnorm;
    % compute the error from each singular triplet
    E = data*V - U*diag(S);
    err_vec = zeros(k,1);
    for j = 1:k
        err_vec(j) = norm(E(:,j));
    end
    err_vecmax = max(err_vec)/data_2norm;
    % output for analysis
    fprintf('%5d %e %e %e \n',k,te,err_mat,err_vecmax)
    clear U S V ts te err_mat E err_vec err_vecmax
end


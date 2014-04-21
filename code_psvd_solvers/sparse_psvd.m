function [U, D, V] = sparse_psvd(A, k, opts)
%function [U, D, V] = sparse_psvd(A, k, opts)
%   This function computes the largest/smallest k singular triplets
%   of matrix A. The U and V are left and right singular vectors of
%   A, and D is the diagonal matrix that stores the singular
%   values. The input options have the default values:
%
%   opts.maxiter    2000
%   opts.tol        1e-5
%   opts.t_start     2*k
%   opts.t_step       20
%   opts.sing       'LA' --- largest singular components
%       another legal value is 'SA' --- smallest
%
%   Note: The option maxiter in some sense heavily depends on the
%   value of k. You might need to tune this parameter such that the
%   maximum number of iterations is not easily exceeded while
%   avoiding allocating more memory than necessary.
%
%   The computation of the SVD is based on the Lanczos algorithm.
%
%   Example:
%      A = sprand(1e5,2e4,1e-4);
%      tic; [U1 D1 V1] = svdk(A, 10); toc
%      tic; [U2 D2 V2] = svds(A, 10); toc
%
%   See also SVDS.

% Copyright 2008, Jie Chen.
% This program is free software; you can redistribute and/or
% modify it for NON-COMMERCIAL purposes. This program is
% distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY, including that of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE.
% $Date: 2008/12/10 18:55:02$

% Modified 04/2013, Zheng Wang

if ~exist('opts','var')      opts = struct;        end
if ~isfield(opts,'maxiter')  opts.maxiter = 2000;  end
if ~isfield(opts,'tol')      opts.tol     = 1e-5;  end
if ~isfield(opts,'t_start')  opts.t_start = 2*k;   end
if ~isfield(opts,'t_step')   opts.t_step  = 20;    end
if ~isfield(opts,'sing')     opts.sing    = 'LA';  end

if strcmp(opts.sing,'LA')==0 && strcmp(opts.sing,'SA')==0
  error('Invalid opts.sing!');
end

[m,n] = size(A);

if m >= n
  
  k = min(k,n);
  opts.maxiter = min(opts.maxiter,n);
  if opts.t_start >= n
    opts.t_start = k+round(k/20);
  end
  
  [V, vecDD] = lanczos(A, k, opts);
  
  vecD = sqrt(vecDD);
  D = diag(vecD);
  U = zeros(m,k);
  for i = 1:k
    U(:,i) = A*V(:,i)/vecD(i);
  end
  
else
  
  k = min(k,m);
  opts.maxiter = min(opts.maxiter,m);
  if opts.t_start >= m
    opts.t_start = k+round(k/20);
  end
  
  [U, vecDD] = lanczos(A', k, opts);
  
  vecD = sqrt(vecDD);
  D = diag(vecD);
  V = zeros(n,k);
  for i = 1:k
    V(:,i) = (U(:,i)'*A)'/vecD(i);
  end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V, vecD] = lanczos(A, k, opts)
% This function returns the largest k eigenpairs of A'A.
% This function is meant to be called by svdk() only.
% The returns, V and vecD, are the eigenvectors and values of A'A.

maxiter = opts.maxiter;
tol     = opts.tol;
t_start = opts.t_start;
t_step  = opts.t_step;
sing    = opts.sing;

n       = size(A,2);

beta    = zeros(maxiter+1,1);
alpha   = zeros(maxiter,1);
Q       = zeros(n,maxiter+1);
qi_1    = zeros(n,1);
T       = zeros(maxiter+1);

Q(:,1)  = randn(n,1);
Q(:,1)  = Q(:,1)/norm(Q(:,1));

i = 1;
% sum_eig_val_prev = -1;
while i <= maxiter
  %--lanczos step
  qi = Q(:,i);
  w = ((A*qi)'*A)' - beta(i)*qi_1;
  alpha(i) = w'*qi;
  w = w - alpha(i)*qi;
  %%----------------------
  %%  re-orthogonalization
  w = w - Q(:,1:i-1)*(w'*Q(:,1:i-1))';
  %%----------------------
  beta(i+1) = norm(w);
  Q(:,i+1) = w/beta(i+1);
  qi_1 = qi;
  iter = i;
  
  T(i,i) = alpha(i);
  T(i,i+1) = beta(i+1);
  T(i+1,i) = beta(i+1);
  
%   %--convergence test
%   if (i>=t_start && mod(i,t_step)==0) || (i==maxiter)
%     eig_val = eig(T(1:iter,1:iter));
%     if strcmp(sing,'LA')
%       [eig_val,idx] = sort(eig_val,'descend');
%     elseif strcmp(sing,'SA')
%       [eig_val,idx] = sort(eig_val,'ascend');
%     end
%     sum_eig_val = sum(eig_val(1:k));
%     if abs(sum_eig_val-sum_eig_val_prev) < tol
%       break;
%     else
%       sum_eig_val_prev = sum_eig_val;
%     end
%   end

  % use a more strict convergence test
  if (i>=t_start && mod(i,t_step)==0) || (i==maxiter)
      [eig_vec, eig_VAL] = eig(T(1:iter,1:iter));
      eig_val = diag(eig_VAL);
      if strcmp(sing,'LA')
          [eig_val,idx] = sort(eig_val,'descend');
          eig_vec = eig_vec(:,idx);
          n_mat = eig_val(1);
      elseif strcmp(sing,'SA')
          [eig_val,idx] = sort(eig_val,'ascend');
          eig_vec = eig_vec(:,idx);
          n_mat = eig_val(iter);
      end
      % check residual norm for convergence
      n_res = abs(beta(iter+1)*eig_vec(iter,1:k+round(k/20))); % residual norm of the first k eigenpairs
      n_conv = sum(n_res < tol*n_mat);
      fprintf('In iteration %5d, %5d eigenpairs converged. \n',iter, n_conv);
      if n_conv == k+round(k/20)
          fprintf('Computation finished in %5d iterations.\n',iter);
          break  % break the while-loop
      end
  end
  i = i + 1;
end

fprintf('maxiter = %d.\n',maxiter)
fprintf('Program stops after iteration %d is finished.\n',i-1)

% %--compute V and vecD
% [VV,DD] = eig(T(1:iter,1:iter));
% if strcmp(sing,'LA')
%   [vecD,idx] = sort(diag(DD),'descend');
% elseif strcmp(sing,'SA')
%   [vecD,idx] = sort(diag(DD),'ascend');
% end
% V = Q(:,1:iter)*VV(:,idx(1:k));
% vecD = vecD(1:k);

% return V and vecD
V = Q(:,1:iter)*eig_vec(:,1:k);
vecD = eig_val(1:k);



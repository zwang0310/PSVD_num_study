function [U S V] = rsvd_k(A, k, p, q, randtype)
% k meand Krylov subspace
%
% k -- the number of singular triplets to compute
% p -- oversampling parameter
% q -- power (basic method if q=0)
% randtype -- distribution underlying the random matrix

% Zheng Wang
% 2014


[m,n] = size(A);

if m >= n   
    F = zeros(m, (k+p)*(1+q));
    % draw a n x (k+p) random testing matrix H
    if strcmp(randtype, 'std_normal')
        H = randn(n,k+p);
    elseif strcmp(randtype, 'uniform')
        H = 2*rand(n,k+p)-ones(n,k+p);
    end
    % apply ((AA')^q)A to the random matrix
    H = A*H; % m x (k+p)
    F(:,1:k+p) = H; 
    for its = 1:q  % power method if q >= 2
        H = (H'*A)';
        H = A*H;
        F(:,(k+p)*its+1:(k+p)*(its+1)) = H;
    end
     
    % compute the orthonormal basis
    [Q,R] = qr(F,0); % Q is m x (k+p)(1+q)

    % solve svd on a reduced matrix Q'A: (k+p)(1+q)xn 
    [U, S, V] = svd(Q'*A, 'econ');
    U = Q*U;

else
    F = zeros(n, (k+p)*(1+q));
    % draw a n x (k+p) random testing matrix H
    if strcmp(randtype, 'std_normal')
        H = randn(k+p,m);
    elseif strcmp(randtype, 'uniform')
        H = 2*rand(k+p,m)-ones(k+p,m);
    end
    % apply ((AA')^q)A to the random matrix
    H = (H*A)'; % n x (k+p)
    F(:,1:k+p) = H;
    for its = 1:q
        H = A*H;
        H = (H'*A)';
        F(:,(k+p)*its+1:(k+p)*(its+1)) = H;
    end
    % compute the orthonormal basis
    [Q,R] = qr(F,0); % Q is n x (k+p)(1+q)

    % solve svd on a reduced matrix AQ: m x(k+p)(1+q)
    [U, S, V] = svd(A*Q, 'econ');
    V = Q*V;    
end


U = U(:,1:k);
V = V(:,1:k);
S = S(1:k,1:k);


end
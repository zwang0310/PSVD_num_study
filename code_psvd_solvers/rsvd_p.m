function [U S V] = rsvd_p(A, k, p, q, randtype)
% p means power method ((A*A')^q)*(A*H)
%
% k -- the number of singular triplets to compute
% p -- oversampling parameter
% q -- power (basic method if q=0)
% randtype -- distribution underlying the random matrix

% Zheng Wang
% 2014


[m,n] = size(A);

if m >= n
    % draw a n x (k+p) random testing matrix H
    if strcmp(randtype, 'std_normal')
        H = randn(n,k+p);
    elseif strcmp(randtype, 'uniform')
        H = 2*rand(n,k+p)-ones(n,k+p);
    end
    % apply ((AA')^q)A to the random matrix
    H = A*H; % m x (k+p)
    for its = 1:q  % power method if q >= 1
        H = (H'*A)';
        H = A*H;
    end
     
    % compute the orthonormal basis
    [Q,R] = qr(H,0); % Q is m x (k+p)

    % solve svd on a reduced matrix Q'A: (k+p)xn 
    [U, S, V] = svd(Q'*A, 'econ');
    U = Q*U;

else
    % draw a n x (k+p) random testing matrix H
    if strcmp(randtype, 'std_normal')
        H = randn(k+p,m);
    elseif strcmp(randtype, 'uniform')
        H = 2*rand(k+p,m)-ones(k+p,m);
    end
    % apply ((AA')^q)A to the random matrix (q = 0)
    H = (H*A)'; % n x (k+p)
    for its = 1:q
        H = A*H;
        H = (H'*A)';
    end
    % compute the orthonormal basis
    [Q,R] = qr(H,0); % Q is n x (k+p)

    % solve svd on a reduced matrix AQ: m x(k+p)
    [U, S, V] = svd(A*Q, 'econ');
    V = Q*V;    
end


U = U(:,1:k);
V = V(:,1:k);
S = S(1:k,1:k);


end
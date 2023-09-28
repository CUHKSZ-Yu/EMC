function [D_new] = correct_emc(D0, epsilon)
% function [D_new] = correct_emc(D0, epsilon)
%
% Correct a distance matrix by Embedding-based Matrix Correction (EMD) method
% via solving a convex optimization problem (ref. Section 3.2)
%
% @param  D0        Initial distance matrix
% @param  epsilon   Default 0.02
% @return D_new     Corrected distance matrix
%
% <Reference>
% [1] Fangchen Yu, et al. "Highly-Efficient Robinson-Foulds Distance Estimation with 
% Matrix Correction", ECAI, 2023.
% [2] Wenye Li, Fangchen Yu. "Calibrating Distance Metrics Under Uncertainty." ECML, 2022.
% [3] Wenye Li. "Scalable calibration of affinity matrices from incomplete
% observations.", ACML, 2020.

if (nargin < 2)
    epsilon = 0.02;
end

r = 10;
gamma = -epsilon / max(D0(:));
low = exp(-epsilon);
S0 = exp(gamma.*D0);
maxiter = 100;
n = size(D0, 1);
np = cell(r, 1);

for t = 1 : maxiter
    p = reshape(mod(randperm(2*n),n)+1, 2*n/r, r);
    for i = 1 : r
        np{i} = unique(p(:,i));
    end
    St = cyclic_cal(S0, np, 2);
    S0 = St;
end
D_new = log(St) ./ gamma;
D_new(1:n+1:n^2) = 0;
D_new = real(D_new);
D_new = (D_new + D_new') / 2;

end


%%
function [Jt] = cyclic_cal(J0, p, iter)
% function [Jt] = cyclic_cal(J0, p, iter)
% 
% @param S0     Initial matrix to calibrate
% @param p      A cell arracy; p{i}: indices of elements in Ci 
% @param iter   # of iterations (default 2)

if (nargin < 3)
    iter = 2;
end
r = length(p);
nc = zeros(r, 1);
It = cell(r,1);
for i = 1 : r
    nc(i) = length(p{i});
    It{i} = zeros(nc(i));
end

Jt = J0;
for t = 1 : iter 
    for i = 1 : r
            P = Jt(p{i}, p{i}) - It{i};
            Q = nearpsd(P, 15);
            Jt(p{i},p{i}) = Q;
            It{i} = Q - P;
    end
end

end


%%
function [X, iter] = nearpsd(A, maxits, low, high)
% function [X, iter] = nearpsd(A, maxits, low, high)
%
% Computes the nearest positive semi-definite matrix 
% for a given square matrix.
%
% @param A        a square matrix to be calibrated
% @param maxits   max num of iters allowed, default 100
% @param low      default 0 
% @param high     default 1
%
% @return X       nearest psd matrix to A
% @return iter    number of iterations taken

if  ~isequal(A,A')
    A = (A + A') / 2;
end
if nargin < 4
    high = 1;
end
if nargin < 3
    low = 0; 
end
if nargin < 2
    maxits = 100;
end

% threshold for convergence & eigs
tolconv = 1.0e-5;
toleigs = 1.0e-5;

n = size(A,1);

U = zeros(n);
Y = A;

[V, D] = eig(Y);
d = diag(D);

iter = 0;
while 1
    T = Y - U;
    
    % project onto psd matrices
    [Q, D] = eig(T);
    d = diag(D);
    p = d > toleigs*d(n);
    X = Q(:,p) * D(p,p) * Q(:,p)';

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxits
        %fprintf('Max iterations reached. ');
        break; 
    end
    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
    
    % problem-dependent here
    Y = X;
    Y(1:n+1:n*n) = 1;
    Y(Y<low) = low;
    Y(Y>high) = high;
end

Y(1:n+1:n*n) = 1;
Y(Y<low) = low;
Y(Y>high) = high;

end
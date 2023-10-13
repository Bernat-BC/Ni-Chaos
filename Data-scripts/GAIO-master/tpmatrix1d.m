% This is a slight adaptation of the original file tpmatrix.m that can be
% found in GAIO-master/matlab. The method and output have not been changed. 


function T_total = tpmatrix1d(tree, f, X, depth, verbose, noise, eps)
 
% TPMATRIX   transition probability matrix.
%
%   TPMATRIX(t, f, X, d, v, n, e) computes the matrix T of transition
%   probabilities between the boxes in the tree t
%
%   t       tree containing the box covering 
%   f       map
%   X       m x d-matrix of sample points
%   d       depth of the tree on which the matrix is computed
%   v       verbose flag: '0' or '1' ('1' prints progress)
%   n       number of noise realisations
%   eps     value of the epsilon inflation. 


d = tree.dim;
b = tree.boxes(depth);                                                     % get the geometry of the boxes
N = size(b,2); S = whos('X'); l = floor(5e7/S.bytes);
T_total = sparse([], [],[], N,N);

if verbose 
    h = waitbar(0,'Computing Transition Probability Matrix...');
end

NOISE = linspace(-eps, eps, noise);
% NOISE = normrnd(0,sqrt(eps), noise);
for k2 = 1:noise
    omega = NOISE(k2);
    IJS = [];
    for k = 0:floor(N/l)                                                   % split in chunks of size l
        K = k*l+1:min((k+1)*l,N);
        c = b(1:d,K); r = b(d+1:2*d,1);                                    % center and radii of the boxes
        n = size(c,2); E = ones(n,1);                                      % sample points in all boxes
        P = kron(E,X)*diag(r) + kron(c',ones(size(X,1),1));
        image = f(P, omega.*ones(1,length(P)));
        I = tree.search(image(:,1)', depth);                               % get box numbers of images
        J = kron(K',ones(size(X,1),1));                                    % column numbers
        pI = find(I>0);                                                    % I <= 0 iff f(P) is outside of covering
        [I,J,S] = find(sparse(I(pI), J(pI), 1, N, N));                     % I = row number >0
                                                                           % J = column number nonzero elements 
                                                                           % S = entries of nonzero elements
        IJS = sparse([IJS; I,J,S]);
    end                                                                    % transition matrix for one noise
    T = sparse(IJS(:,1), IJS(:,2), IJS(:,3), N, N);   
    T_total = T_total + T;                         
    if verbose, waitbar(k2/noise); end
end
if verbose
    try
        close(h); 
    catch
        h = waitbar(0,'Waitbar');
        close(h);
    end
end
cs = ones(2^depth);
T_total = sparse(T_total*spdiags(1./cs', 0, N, N)./(noise*length(X)));     % normalize
if verbose, fprintf('\n'); end
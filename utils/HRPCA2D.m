function [L,S] = HRPCA2D(Y,varargin)
% =========================
%  haart2 + 4-nn
% =========================
% =========================
% parse the input arguments
% =========================
sizeY = size(Y);
sizeD = [sizeY(1)*sizeY(2), sizeY(3)];
D = reshape(Y, sizeD); %clear Y

[maxIter, tol, rho, lambda, mu_a, mu_b, gamma, debug, S, L, print_freq, GT] = ...
    process_options(varargin, ...
    'maxIter',   100, ...
    'tol',       1e-3, ...
    'rho',       1.05, ...
    'lambda',    2, ...
    'mu_a',      4, ...
    'mu_b',     2, ...
    'gamma',     [1,1,1,1], ...
    'debug',     1, ...
    'S',         zeros(sizeY), ...
    'L',         Y, ...
    'print_freq',5, ...
    'GT',        NaN ...
    );

% set hyper-parameter for sparse penalty
lambda = lambda/sqrt(sizeY(1)*sizeY(2));

% Initialize mu
normD   = norm(D,'fro');
[~,norm_two,~] = svdsecon(D, 1);
norm_inf = max(abs(D(:)))/lambda;
dual_norm = max(norm_two, norm_inf);
mu_a = mu_a/dual_norm; % this one can be tuned
max_mu_a = mu_a * 1e8;
mu_b = mu_b/dual_norm;
max_mu_b = mu_b * 1e8;

% Initialize multiplier
A  = zeros(sizeY);
B1 = 0;
B2 = 0;
B3 = 0;
B4 = 0;

% Iteration
for iter = 1:maxIter
    % Update S (Sparse noise)
    S = softthre(Y-L+A/mu_a, lambda/mu_a); % 3D data

    % Update C (Haar coefficients)
    [FL1, FL2, FL3, FL4] = haar_trans2(L);  % forward Haar Transform  % 3D data
    sizeF = size(FL1);

    FL1 = reshape(FL1, sizeF(1)*sizeF(2), sizeF(3));
    FL2 = reshape(FL2, sizeF(1)*sizeF(2), sizeF(3));
    FL3 = reshape(FL3, sizeF(1)*sizeF(2), sizeF(3));
    FL4 = reshape(FL4, sizeF(1)*sizeF(2), sizeF(3)); % 2D data

    C1 = prox_nn(FL1+B1/mu_b, gamma(1)/mu_b);
    C2 = prox_nn(FL2+B2/mu_b, gamma(2)/mu_b);
    C3 = prox_nn(FL3+B3/mu_b, gamma(3)/mu_b);
    C4 = prox_nn(FL4+B4/mu_b, gamma(4)/mu_b); % 2D data

    CB1 = reshape(C1-B1/mu_b, sizeF);
    CB2 = reshape(C2-B2/mu_b, sizeF);
    CB3 = reshape(C3-B3/mu_b, sizeF);
    CB4 = reshape(C4-B4/mu_b, sizeF); % 3D data

    % Update L
    rhs = mu_a*(Y-S+A/mu_a);
    rhs = rhs + mu_b*ihaar_trans2(CB1,CB2,CB3,CB4,sizeY(1),sizeY(2));
    L = rhs/(mu_a+mu_b);

    %myshow(L)

    % stop criterion
    stopC = norm(reshape(Y-L-S, sizeD),'fro')/normD;
    converge_flag = stopC<tol;
    if isnan(GT)
        if debug && (mod(iter,print_freq)==0 || iter==1 || iter==maxIter || converge_flag)
            disp(['iter ' num2str(iter) ',mu=' num2str(mu_a,'%2.1e')  ...
                ',Y-X-S=' num2str(stopC,'%2.3e')]);
        end
    else
        if debug && (mod(iter,print_freq)==0 || iter==1 || iter==maxIter || converge_flag)
            disp(['iter ' num2str(iter) ',mu=' num2str(mu_a,'%2.1e')  ...
                ',Y-X-S=' num2str(stopC,'%2.3e'),',PSNR=' num2str(psnr(GT,L),'%2.2f')]);
        end
    end

    % Update multiplier
    if converge_flag
        break;
    else
        A = A + mu_a*(Y-L-S);
        B1 = B1 + mu_b*(FL1-C1);
        B2 = B2 + mu_b*(FL2-C2);
        B3 = B3 + mu_b*(FL3-C3);
        B4 = B4 + mu_b*(FL4-C4);
        mu_a = min(max_mu_a,mu_a*rho);
        mu_b = min(max_mu_b,mu_b*rho);
    end

end


end
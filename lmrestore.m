function [tn, flag] = lmrestore(lmobj, b)
% function [tn, flag] = lmrestore(lmobj, b)
% regularized image restoration for lm 
% Input: 
%     lmobj - lm data. see lmdatainit
%     b - weight factor for the regularization term
% Output:
%     tn - restored image
%     flag - 0 normal exit. 1 - too many iterations, 2-line search stepsize
%     too small

maxiter = 1000;
maxiter_l = 100;
iscale=1e2;
flag = 0;
[w,h]=size(lmobj.initimg);
t = 0.5;

stepsize_limit = 1e-4 / ((h*w)^2);

b1 = [1 -2 1]';
b2 = [1 -2 1];
b3 = sqrt(2) .* [1 -1 0;-1 1 0; 0 0 0];
b1r = b1;
b2r=b2;
b3r = sqrt(2) .* [0 0 0; 0 1 -1; 0 -1 1];
g=zeros(w,h,3);

theta = ones(w,h);
theta = lmdeconv(lmobj, 50, theta);
theta = theta / sum(theta(:));

it = 1;

while (it < maxiter)
    U = lmdeconv(lmobj,1,theta);
    % compute gradient
    g(:,:,1) = filter2(b1, theta);
    g(:,:,2) = filter2(b2, theta);
    g(:,:,3) = filter2(b3, theta);
    g2 = sum(g.^2,3) + theta.^2/iscale;
    g = g ./ (g2 + eps);
    g(:,:,1) = filter2(b1r, g(:,:,1));
    g(:,:,2) = filter2(b2r, g(:,:,2));
    g(:,:,3) = filter2(b3r, g(:,:,3));
    g0 = b .* theta .* (sum(g,3) + theta / iscale ./ (g2 + eps) - h*w) - U + sum(U(:)).*theta;
    %g0 = b .* theta .* (sum(g,3) + theta./(g2 + eps) - nnz(g2) + sum(U(:))/b) - U;
    f0 = 0.5 * b * sum(log(g2(:)+eps)) - sum(U(:).*log(theta(:)+eps));

    % line search 
    a0 = 1/sum(U(:));
    f1 = f0 + 1;
    it_l = 1;
    while (f1 > f0 && it_l < maxiter_l)
        tn = theta .* exp(-g0*a0);
        tn = tn / sum(tn(:));
        g(:,:,1) = filter2(b1, tn);
        g(:,:,2) = filter2(b2, tn);
        g(:,:,3) = filter2(b3, tn);
        g2 = sum(g.^2,3) + tn.^2 / iscale;
        f1 = 0.5 * b * sum(log(g2(:)+eps)) - sum(U(:) .* log(tn(:) +eps));
        a0 = a0 * t;
        it_l = it_l + 1;
    end
    
    if (it_l >= maxiter_l)
        flag = 2;
        break;
    end

    if (norm(tn(:) - theta(:)) < stepsize_limit)
        flag = 0;
        break;
    end

    theta = tn;
    it = it + 1;
    
    if (mod(it,10)==0)
        disp(['iter: ' int2str(it)]);
    end
end

if (it >= maxiter)
    flag = flag + 1;
end

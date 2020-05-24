function img = lmdeconv(lmobj, iters, prev, prior)
% img = lmdeconv(locdata, iters, prev, prior)
% Perform statitically accurate rendering of SMLM data using ML estimator
% 
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iter - # of iterations
%   prev - result from previos call of this function. Optional
%   priot - Dirichlet prior, Optional.
% Outputs -
%   img - deconv img

narginchk(2,4);

if(~isscalar(iters))
    error('iters must be a scalar');
end
if(~exist('prior','var'))
    prior=1;
end
if(~exist('prev','var') || isempty(prev))
    prev = ones(lmobj.imgsize);
end

for i = 1:iters
    if (length(lmobj.imgsize) == 2)
        prev = lmdeconvmex(lmobj.data, prev, lmobj.psfs, double(prior));
    else
        prev = lmdeconv3dmex(lmobj.data, prev, lmobj.psfs, lmobj.zpsfs, double(prior));
    end
    if (mod(i,10) == 0)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

img = prev;

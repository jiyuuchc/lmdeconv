function img = lmdeconv(lmobj, iters, prev, prior)
% img = lmdeconv(locdata, iters, prev)
% perform deconvolution on SMLM data
%
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iter - # of iterations
%   prev - result from previos iteration
%
% Outputs -
%   img - deconv img

narginchk(2,4);

if(~isscalar(iters))
    error('iters must be a scalar');
end

if(~exist('prev','var') || isempty(prev))
    prev = lmobj.initimg;
end

if(~exist('prior','var'))
    prior=1;
end

for i = 1:iters
    prev = lmdeconvmex(lmobj.data, prev, lmobj.psfs, double(prior));
    if (mod(i,10) == 0)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

img = prev;

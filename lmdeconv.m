function img = lmdeconv(lmobj, iters, prev)
% img = lmdeconv(locdata, iters, prevsample)
% perform deconvolution on SMLM data
%
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iter - # of iterations
%   prev - result from previos iteration
%
% Outputs -
%   img - deconv img

narginchk(2,3);

if(~isscalar(iters))
    error('iters must be a scalar');
end

if(~exist('pre','var') || isempty(prev))
    prev = lmobj.img;
end

for i = 1:iters
    prev = lmdeconvmex(lmobj.data, prev, lmobj.psfs);
    if (mod(i,10) == 0)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

img = prev;

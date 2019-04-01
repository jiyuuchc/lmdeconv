function samples = lmsample(lmobj, iters, prev, skip, prior)
% samples = lmdsample(locdata, iters, prev, skip)
% Sample the posterior distribution on theta based on SMLM data
%
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iter - # of samples to be drawn
%   prev - Optional. result from previos iteration. 
%   skip - Optional. skipping samples, default 1
%   prior - Dir prior. can be a scaler or a matrix.
%
% Outputs -
%   samples - samples

narginchk(2,5);

if (~exist('prior','var'))
    prior = 1.0;
end

if (~exist('skip','var'))
    skip = 1;
else
    if (~isscalar(skip))
        error ('skip must be a scalar');
    end
end

if(~isscalar(iters))
    error('iters must be a scalar');
end

if(~exist('prev','var') || isempty(prev))
    prev = lmobj.initimg;
end

if (~isscalar(prior) && numel(prior) ~= numel(prev))
    error('Prior size does not match image');
end

[h,w] = size(prev);
samples = zeros(h,w,iters);
for i = 1:iters * skip
    prev = lmsamplermex(lmobj.data, prev, lmobj.psfs, double(prior));
    if (mod(i,skip)==0)
        samples(:,:,i/skip) = prev;
    end
    if (mod(i,100) == 0)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

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
    prev = histimg(lmobj);
end

if (~exist('prior','var'))
    prior = 1/2; % Jeffreys prior
end
if (~isscalar(prior) && numel(prior) ~= numel(prev))
    error('Prior size does not match image');
end


samplesize = [size(prev) iters];
samples = zeros(samplesize);
for i = 1:iters * skip
    if (length(lmobj.imgsize) == 2)
        prev = lmsamplermex(lmobj.data, prev, lmobj.psfs, double(prior));
        if (mod(i,skip)==0)
            samples(:,:,i/skip) = prev;
        end
    else
        prev = lmsampler3dmex(lmobj.data, prev, lmobj.psfs, lmobj.zpsfs, double(prior));
        if (mod(i,skip)==0)
            samples(:,:,:,i/skip) = prev;
        end
    end
    if (mod(i,1000) == 0)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

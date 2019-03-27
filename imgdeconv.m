function img = imgdeconv(data, psf, iter, scale, prior, prev)
% img = imgdeconv(data, psf, scale, iter, prior)
% perform image deconvolution 
%
% Inputs -
%   data - image data
%   psf - psf
%   iter - # of iterations
%   scale - pixel subsampling 
%   prev - result from previos iteration
%
% Outputs -
%   img - deconv img

narginchk(3,6);
if (~exist('scale','var') || isempty(scale))
    scale = 1;
end

data = uint32(imresize(data,scale)/scale/scale);

if (~exist('prev','var') || isempty(prev))
    prev = double(data);
end

if (exist('prior', 'var') && ~isempty(prior))
    if (~isscalar(prior) && numel(prior) ~= numel(data))
        error('Wrong prior size');
    end
else
    prior = 1;
end

psf = double(psf);
[h,w] = size(data);
[ph,pw] = size(psf);
if (ph ~= pw || mod(ph,2)~=1)
    error('psf must be odd size square matrix');
end

if (h ~= size(prev,1) || w ~= size(prev,2))
    error('wrong prev size');
end

for i = 1:iter
    prev = lmdeconvmex(data, prev, psf, prior);
end
img = prev;

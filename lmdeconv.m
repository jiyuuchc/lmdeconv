function img = lmdeconv(locdata, iters, pixelsize, prevsample)
% img = lmdeconv(locdata, iters, pixelsize, prevsample)
% perform deconvolution on SMLM data
%
% Inputs -
%   locdata - 3 x N matrix. X, Y, Sigma. In same unit (e.g. nm)
%   iter - # of iterations
%   pixlesize - the pixel size of the image to be rendered
%   prevsample - a sample from previous iterations. Default to a histogram
%                of locdata
%
% Outputs -
%   img - deconv img

narginchk(2,4);

numsigmabins = 50;

if (~exist('pixelsize','var'))
    pixelsize = 8;
end

if(~isscalar(iters))
    error('iters must be a scalar');
end

if(~isscalar(pixelsize))
    error('pixelsize must be a scalar');
end

[sigmabins, sigmaedges] = discretize(locdata(3,:),numsigmabins);
psfhsize = max(ceil(sigmaedges(end) * 5 / pixelsize), 10);
psfsize = psfhsize * 2 + 1;
psfs = zeros(psfsize, psfsize, length(sigmaedges)-1);
for i = 1:size(psfs,3)
    psfs(:,:,i) = fspecial('gaussian', psfsize, 0.5 /pixelsize * (sigmaedges(i) + sigmaedges(i+1)));
end

xedges = floor(min(locdata(1,:))/pixelsize)*pixelsize:pixelsize:max(locdata(1,:))+pixelsize;
yedges = floor(min(locdata(2,:))/pixelsize)*pixelsize:pixelsize:max(locdata(2,:))+pixelsize;
xbins = discretize(locdata(1,:),xedges);
ybins = discretize(locdata(2,:),yedges);

if(~exist('prevsample','var') || isempty(prevsample))
    img = histcounts2(locdata(2,:), locdata(1,:), yedges, xedges);
    img = padarray(img, [psfhsize, psfhsize], 'both');
else
    img = prevsample;
end

data = zeros(3, size(locdata,2), 'uint32');
data(1,:) = xbins + psfhsize - 1;
data(2,:) = ybins + psfhsize - 1;
data(3,:) = sigmabins - 1;

for i = 1:iters
    img = lmdeconvmex(uint32(data), double(img), double(psfs));
    if (mod(i,10) == 0)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

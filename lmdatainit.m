function dataobj = lmdatainit(locdata, pixelsize, edges, numsigmabins)
% dataobj = lmdatainit(locdata, pixelsize, edges, numsigmabins)
% Prepare data to be used for other fuctions
%
% Inputs -
%   locdata - 3 x N matrix. X, Y, Sigma. In same unit (e.g. nm)
%   pixlesize - the pixel size of the image to be rendered
%   edges - Optional. Bin edges used to discretize localization data
%   numsigmabins - Optional. default 50.
%
% Outputs -
%   img - deconv img

narginchk(1,4);

if (~exist('numsigmabins','var'))
    numsigmabins = 50;
end

if(~isscalar(numsigmabins))
    error('numsigmabins must be a scalar');
end

if (~exist('pixelsize','var'))
    pixelsize = 8;
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

if (~exist('edges','var'))
    xedges = floor(min(locdata(1,:))/pixelsize)*pixelsize:pixelsize:max(locdata(1,:))+pixelsize;
    yedges = floor(min(locdata(2,:))/pixelsize)*pixelsize:pixelsize:max(locdata(2,:))+pixelsize;
else
    xedges = edges{1} * pixelsize;
    yedges = edges{2} * pixelsize;
end
xbins = discretize(locdata(1,:),xedges);
ybins = discretize(locdata(2,:),yedges);

img = histcounts2(locdata(2,:), locdata(1,:), yedges, xedges);
img = double(padarray(img, [psfhsize, psfhsize], 'both'));

data = zeros(3, size(locdata,2), 'uint32');
data(1,:) = xbins + psfhsize - 1;
data(2,:) = ybins + psfhsize - 1;
data(3,:) = sigmabins - 1;

dataobj = struct();
dataobj.data = uint32(data);
dataobj.initimg = double(img);
dataobj.psfs = double(psfs);
dataobj.pixelsize = pixelsize;

function dataobj = lmdatainit(varargin)
% dataobj = lmdatainit(locdata, pixelsize, padding, numsigmabins)
% dataobj = lmdatainit(locdata, edges, padding, numsigmabins)
% Prepare data to be used for other fuctions
%
% Inputs -
%   locdata - 3 x N matrix. X, Y, Sigma. In same unit (e.g. nm)
%   pixlesize - the pixel size of the image to be rendered
%   edges - Bin edges used to discretize localization data.
%   padding - Optional. default is using the largest psf size
%   numsigmabins - Optional. default 50.
%
% Outputs -
%  dataobj - Object encapsulating necessary data for further analyses

narginchk(2,4);
locdata = double(varargin{1});
if (size(locdata,1) < 3)
    error('input data should have 3 rows');
end

if (length(varargin) < 4)
    numsigmabins = 50;
else
    numsigmabins = varargin{4};
end

if(~isscalar(numsigmabins))
    error('numsigmabins must be a scalar');
end

if (length(varargin) >= 3)
    padding = varargin{3};
    if (isscalar(padding))
        padding = [padding padding];
    end
    if (length(padding) ~= 2)
        error('Wrong padding size');
    end
end
 
if (isscalar(varargin{2}))
    pixelsize = varargin{2};
else
    edges = varargin{2};
    if (~iscell(edges))
        error('edges should be a 2-element cell array');
    end
    xedges = edges{1};
    yedges = edges{2};
    dxedges = diff(xedges);
    dyedges = diff(yedges);
    if (range(dxedges) >1e-10  || range(dyedges) >1e-10 || dxedges(1) - dyedges(1)>1e-10)
        error('edges should be equal spacing');
    end
    pixelsize = double(dxedges(1));
end

[sigmabins, sigmaedges] = discretize(double(locdata(3,:)),numsigmabins);
psfhsize = max(ceil(sigmaedges(25) * 3 / pixelsize), 10);
psfsize = psfhsize * 2 + 1;
psfs = zeros(psfsize, psfsize, length(sigmaedges)-1);
for i = 1:size(psfs,3)
    psfs(:,:,i) = fspecial('gaussian', psfsize, 0.5 /pixelsize * (sigmaedges(i) + sigmaedges(i+1)));
end

if (~exist('xedges','var'))
    xedges = floor(min(locdata(1,:))/pixelsize)*pixelsize:pixelsize:max(locdata(1,:))+pixelsize;
    yedges = floor(min(locdata(2,:))/pixelsize)*pixelsize:pixelsize:max(locdata(2,:))+pixelsize;
end
xbins = discretize(locdata(1,:),xedges);
ybins = discretize(locdata(2,:),yedges);

img = histcounts2(locdata(2,:), locdata(1,:), yedges, xedges);

if (~exist('padding','var'))
    padding = [psfhsize, psfhsize];
end
img = double(padarray(img, padding, 'both'));

data = zeros(3, size(locdata,2), 'uint32');
data(1,:) = xbins + padding(1) - 1;
data(2,:) = ybins + padding(2) - 1;
data(3,:) = sigmabins - 1;

dataobj = struct();
dataobj.data = uint32(data);
dataobj.initimg = double(img);
dataobj.psfs = double(psfs);
dataobj.pixelsize = pixelsize;
dataobj.origin = locdata;
dataobj.edges = {xedges, yedges};
dataobj.padding=padding;

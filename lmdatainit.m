function dataobj = lmdatainit(varargin)
% dataobj = lmdatainit(locdata, pixelsize, padding, numsigmabins)
% dataobj = lmdatainit(locdata, edges, padding, numsigmabins)
% Prepare data to be used for other fuctions
%
% Inputs -
%   locdata - 3 x N matrix. X, Y, Sigma. In same unit (e.g. nm)
%   pixlesize - the pixel size of the image to be rendered
%   edges - Bin edges used to discretize localization data.
%   padding - Optional. 
%   numsigmabins - Optional. default 50.
%
% Outputs -
%  dataobj - Object encapsulating necessary data for further analyses

narginchk(2,4);
locdata = double(varargin{1});
if (size(locdata,1) < 3)
    error('input data should have 3 rows');
end

if (size(locdata,1) == 5) 
    is3d = true;
else
    is3d = false;
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
 
if (~iscell(varargin{2}))
    tmp = varargin{2};
    pixelsize = tmp(1);
    xedges = floor(min(locdata(1,:))/pixelsize)*pixelsize:pixelsize:max(locdata(1,:))+pixelsize;
    yedges = floor(min(locdata(2,:))/pixelsize)*pixelsize:pixelsize:max(locdata(2,:))+pixelsize;    
    if (is3d)
        zpixelsize = tmp(2);
        zedges = floor(min(locdata(4,:))/zpixelsize)*zpixelsize:zpixelsize:max(locdata(4,:))+zpixelsize;
    end
else
    edges = varargin{2};
    if (~iscell(edges))
        error('edges should be a cell array');
    end
    xedges = edges{1};
    yedges = edges{2};
    dxedges = diff(xedges);
    dyedges = diff(yedges);
    if (range(dxedges) >1e-10  || range(dyedges) >1e-10 || dxedges(1) - dyedges(1)>1e-10)
        error('edges should be equal spacing');
    end
    pixelsize = double(dxedges(1));
    if (is3d)
        zedges = edges{3};
        dzedges = diff(zedges);
        if (range(dzedges) > 1e-10)
            error('edges should be equal spacing');
        end
        zpixelsize = double(dzedges(1));
    end
end

xbins = discretize(locdata(1,:),xedges);
ybins = discretize(locdata(2,:),yedges);
if (is3d)
    zbins = discretize(locdata(4,:), zedges);
end

[sigmabins, sigmaedges] = discretize(double(locdata(3,:)),numsigmabins);
psfhsize = max(ceil(sigmaedges(end) * 2 / pixelsize), 10);
psfsize = psfhsize * 2 + 1;
psfs = zeros(psfsize, psfsize, length(sigmaedges)-1);
for i = 1:size(psfs,3)
    psfs(:,:,i) = fspecial('gaussian', psfsize, 0.5 /pixelsize * (sigmaedges(i) + sigmaedges(i+1)));
end

if (is3d)
    [zsigmabins, zsigmaedges] = discretize(double(locdata(5,:)),numsigmabins);
    psfzhsize = max(ceil(zsigmaedges(end) * 2 / zpixelsize), 10);
    zpsfsize = psfzhsize * 2 + 1;
    zpsfs = zeros(zpsfsize, length(sigmaedges)-1);
    for i = 1:size(psfs,3)
        zpsfs(:,i) = sum(fspecial('gaussian', zpsfsize, 0.5 /zpixelsize * (zsigmaedges(i) + zsigmaedges(i+1))),2);
        %zpsfs(:,i) = normpdf((-psfzhsize:psfzhsize) * zpixelsize, 0, (zsigmaedges(i)+zsigmaedges(i+1))/2);
    end
end

if (~exist('padding','var'))
    if (~is3d)
        padding = [psfhsize, psfhsize];
    else
        padding = [psfhsize, psfhsize, psfzhsize];
    end
end

%data = zeros(3, size(locdata,2), 'uint32');
data(1,:) = xbins + padding(1) - 1;
data(2,:) = ybins + padding(2) - 1;
data(3,:) = sigmabins - 1;
if (is3d)
    data(4,:) = zbins + padding(3) -1;
    data(5,:) = zsigmabins - 1;
end

%img = histcounts2(locdata(2,:), locdata(1,:), yedges, xedges);
%img = double(padarray(img, padding, 'both'));

dataobj = struct();
dataobj.data = uint32(data);
%dataobj.initimg = double(img);
dataobj.psfs = double(psfs);
dataobj.pixelsize = pixelsize;
dataobj.origin = locdata;
dataobj.edges = {xedges, yedges};
dataobj.padding=padding;
dataobj.imgsize = [length(yedges)-1+padding(2)*2, length(xedges)-1+padding(1)*2];

if (is3d)
    dataobj.zpsfs = zpsfs;
    dataobj.pixelsize = [pixelsize, zpixelsize];
    dataobj.edges = {xedges, yedges,zedges};
    dataobj.imgsize = [dataobj.imgsize, length(zedges)-1+padding(3)*2];
end

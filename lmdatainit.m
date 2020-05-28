function dataobj = lmdatainit(locdata, template, varargin)
% dataobj = lmdatainit(locdata, pixelsize ...)
% dataobj = lmdatainit(locdata, template)
% Prepare data to be used for other fuctions
%
% Inputs -
%   locdata - 3 x N matrix. X, Y, Sigma. In same unit (e.g. nm)
%   pixlesize - the pixel size of the image to be rendered
%   template - create new dataobj using the pixel arrangement of and old
%              dataobj template
%   Optional parameters (these are ignored in the 2nd calling form):
%       padding - Optional. default to half of the biggest psf. 
%       numsigmabins - Optional. default 50.
%
% Outputs -
%   dataobj - Object encapsulating necessary data for further analyses

parser = inputParser;
addRequired(parser, 'data', @(p) ismatrix(p) && (size(p,1)==3 || size(p,2)==3));
addRequired(parser, 'template', @(p) (isnumeric(p) && isscalar(p)) || isstruct(p));
addParameter(parser, 'padding', [], @(p) isscalar(p) || (isvector(p) && length(p)==2));
addParameter(parser, 'nsigmabins', 50, @(p) isscalr(p) && p > 0);

parse(parser, locdata, template, varargin{:});

locdata = double(locdata);
if size(locdata,1) ~= 3
    locdata = locdata';
end

if ~isstruct(template) % 2nd input is pixelsize
    pixelsize = template;
    
    xedges = floor(min(locdata(1,:))/pixelsize)*pixelsize:pixelsize:max(locdata(1,:))+pixelsize;
    yedges = floor(min(locdata(2,:))/pixelsize)*pixelsize:pixelsize:max(locdata(2,:))+pixelsize;
    xbins = discretize(locdata(1,:),xedges);
    ybins = discretize(locdata(2,:),yedges);
    
    [sigmabins, sigmaedges] = discretize(double(locdata(3,:)),parser.Results.nsigmabins);

    psfs = cell(length(sigmaedges)-1,1);
    for i = 1:size(psfs,1)
        sigma = 0.5 /pixelsize * (sigmaedges(i) + sigmaedges(i+1));
        psfhsize = round(sigma * 1.5);
        psfsize = psfhsize * 2 + 1;
        psfs{i} = fspecial('gaussian', psfsize, sigma);
    end
    
    padding = parser.Results.padding;
    if isempty(padding)
        padding = [psfhsize, psfhsize];
    else
        if isscalar(padding)
            padding = [padding padding];
        end
    end

    data = zeros(3, length(xbins));
    data(1,:) = xbins + padding(1) - 1;
    data(2,:) = ybins + padding(2) - 1;
    data(3,:) = sigmabins - 1;
    
    dataobj = struct();
    dataobj.data = uint32(data);
    dataobj.psfs = psfs;
    dataobj.pixelsize = pixelsize;
    dataobj.padding=padding;
    dataobj.imgsize = [length(yedges)-1+padding(2)*2, length(xedges)-1+padding(1)*2];
    dataobj.X = (0:dataobj.imgsize(2)-1)*pixelsize + xedges(1) - padding(1) * pixelsize + pixelsize /2;
    dataobj.Y = (0:dataobj.imgsize(1)-1)*pixelsize + yedges(1) - padding(2) * pixelsize + pixelsize /2;
    dataobj.S = (sigmaedges(1:end-1) + sigmaedges(2:end))/2;
else % 2nd input is a dataobj
    dataobj = template;
    xedges = dataobj.X - dataobj.pixelsize / 2;
    yedges = dataobj.Y - dataobj.pixelsize / 2;

    data = zeros(3, size(locdata,2));
    data(1,:) = discretize(locdata(1,:),xedges) - 1;
    data(2,:) = discretize(locdata(2,:),yedges) - 1;
    if ( any(isnan(data(1,:))) || any(isnan(data(2,:))) ...
            || min(data(1,:)) < dataobj.padding(1) ...
            || max(data(1,:)) >= dataobj.imgsize(2) - dataobj.padding(1) ...
            || min(data(2,:)) < dataobj.padding(2) ...
            || max(data(2,:)) >= dataobj.imgsize(1) - dataobj.padding(2))
        error('XY data out of range');
    end
    
    ds = (dataobj.S(2) - dataobj.S(1))/2;
    sigmaedges = [dataobj.S - ds, dataobj.S(end) + ds];
    data(3,:) = discretize(locdata(3,:), sigmaedges) - 1;
    if (any(isnan(data(3,:))))
        error('Sigma out of range');
    end
    dataobj.data = uint32(data);
end

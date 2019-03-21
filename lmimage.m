function img=lmimage(locdata, pixelsize)

xedges = floor(min(locdata(1,:))/pixelsize)*pixelsize:pixelsize:max(locdata(1,:))+pixelsize;
yedges = floor(min(locdata(2,:))/pixelsize)*pixelsize:pixelsize:max(locdata(2,:))+pixelsize;

psfhsize = max(ceil(max(locdata(3,:)) * 5 / pixelsize), 10);
psfsize = psfhsize * 2 + 1;

xbins = discretize(locdata(1,:),xedges) + psfhsize;
ybins = discretize(locdata(2,:),yedges) + psfhsize;
sigmas = locdata(3,:) / pixelsize;
img = zeros(length(yedges) + psfsize - 2, length(xedges) + psfsize -2);
for i = 1:size(locdata,2)
    img(ybins(i)-psfhsize:ybins(i)+psfhsize,xbins(i)-psfhsize:xbins(i)+psfhsize) = img(ybins(i)-psfhsize:ybins(i)+psfhsize,xbins(i)-psfhsize:xbins(i)+psfhsize) + fspecial('gaussian',psfsize, sigmas(i));
end
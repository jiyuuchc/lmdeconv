function img=gausimg(lmobj)
% img = gausimg(lmobj)
% render Gaussian spot image using lmobj as input
% input:
%   lmobj: SMLM data object (see lmdatainit )
% output:
%   img: Rendered image.

img = zeros(lmobj.imgsize);
[h0,w0]=size(img);
for i = 1:length(lmobj.data)
    x0 = lmobj.data(1,i)+1;
    y0 = lmobj.data(2,i)+1;
    p0 = lmobj.data(3,i)+1;
    psf = lmobj.psfs{p0};
    [h,w]=size(psf);
    wh = (w-1)/2;
    hh = (h-1)/2;
    if (x0 > wh && x0+wh <= w0 && y0 > hh && y0+hh <= h0)
        img(y0-hh:y0+hh,x0-wh:x0+wh) = img(y0-hh:y0+hh,x0-wh:x0+wh) + psf;
    end
end

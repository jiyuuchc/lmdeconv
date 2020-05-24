function img=gausimg(lmobj)
% img = gausimg(lmobj)
% render Gaussian spot image using lmobj as input
% input:
%   lmobj: SMLM data object (see lmdatainit )
% output:
%   img: Rendered image.

img = zeros(size(lmobj.initimg));
[h0,w0]=size(img);
[h,w,~]=size(lmobj.psfs);
wh = (w-1)/2;
hh = (h-1)/2;
for i = 1:length(lmobj.data)
    x0 = lmobj.data(1,i)+1;
    y0 = lmobj.data(2,i)+1;
    p0 = lmobj.data(3,i)+1;
    if (x0 > wh && x0+wh <= w0 && y0 > hh && y0+hh <= h0)
        img(y0-hh:y0+hh,x0-wh:x0+wh) = img(y0-hh:y0+hh,x0-wh:x0+wh) + lmobj.psfs(:,:,p0);
    end
end

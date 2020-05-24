function img = histimg(lmobj)
% img = histimg(lmobj)
% render histogram representation of image using lmobj as input
% input:
%   lmobj: SMLM data object (see lmdatainit )
% output:
%   img: Rendered image.

img = zeros(lmobj.imgsize);
if (length(lmobj.imgsize) == 2)
    for i = 1:size(lmobj.data, 2)
        img(lmobj.data(2,i)+1,lmobj.data(1,i)+1) = img(lmobj.data(2,i)+1,lmobj.data(1,i)+1) + 1;
    end
else
    for i = 1:size(lmobj.data, 2)
        img(lmobj.data(2,i)+1,lmobj.data(1,i)+1, lmobj.data(4,i)+1) = img(lmobj.data(2,i)+1,lmobj.data(1,i)+1,lmobj.data(4,i)+1) + 1;
    end
end

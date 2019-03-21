%load('testdata_cross');
load('testdata_ring');
space = [400, 80, 20];
figure;
colormap('gray');
subplot(4,3,1);
imagesc(img0);
axis equal; axis off; 

for i = 1:3

    data0=data(:,1:space(i):end);
    lmobj = lmdatainit(data0, 3.2);
    disp(['Num of localizations: ' int2str(length(data0))]);
    disp('Start reconstruction...');
    
    tic
    img = lmdeconv(lmobj, 30);
    toc

    subplot(4,3,1 + i*3);
    imagesc(lmimage(data0, 3.2));
    axis equal; axis off;

    subplot(4,3,2 + i*3);
    imagesc(deconvlucy(lmobj.initimg, fspecial('gaussian',30,2.2187)));
    axis equal; axis off;

    subplot(4,3,3 + i*3);
    imagesc(img);
    axis equal;axis off;

end

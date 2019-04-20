%load('testdata_cross');
load('testdata_ring');
space = [ 2000, 600, 200];
figure;
colormap('gray');
subplot(4,3,1);
imagesc(img0(473-139:473+139, 473-139:473+139));
axis equal; axis off; 
edges{1}=-40:1.6:40;
edges{2}=-40:1.6:40;
for i = 1:3

    data0=data(:,1:space(i):end);
    lmobj = lmdatainit(data0, edges, 45);
    disp(['Num of localizations: ' int2str(length(data0))]);
    disp('Start reconstruction...');
    
    tic
    img = lmdeconv(lmobj, 50);
    toc

    subplot(4,3,1 + i*3);
    imagesc(crop(lmobj.initimg));
    axis equal; axis off;

    subplot(4,3,2 + i*3);
    imagesc(crop(deconvlucy(lmobj.initimg, fspecial('gaussian',31,4.46))));
    axis equal; axis off;

    subplot(4,3,3 + i*3);
    imagesc(crop(img));
    axis equal;axis off;

end


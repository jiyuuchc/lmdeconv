clear;
load('testdata_ring');
space = [2000, 600, 200];
figure;
colormap('gray');
subplot(4,3,1);
imagesc(img0(473-139:473+139, 473-139:473+139));
axis equal; axis off; 
edges{1}=-40:1.6:40;
edges{2}=-40:1.6:40;
ntests = length(space);
for i = 1:ntests

    data0=data(:,1:space(i):end);
    disp(['Test ' int2str(i) ' - ' int2str(length(data0)) ' localizations.']);

    lmobj = lmdatainit(data0, edges, 45);

    % plot original data
    subplot(ntests+1,3,1 + i*3);
    imagesc(crop(lmobj.initimg));
    axis equal; axis off;


    % compute and display MAP estimator
    disp('MAP reconstruction...');
    tic
    img = lmdeconv(lmobj, 50);
    toc
    subplot(ntests+1,3,2 + i*3);
    imagesc(crop(img));
    axis equal;axis off;

    
    % compute and display MMSE estimator
    disp('MMSE reconstruction...');
    tic
    samples = lmsample(lmobj, 10000);
    toc
    img = mean(samples(:,:,5000:end),3);
    subplot(ntests+1,3,3 + i*3);
    %imagesc(crop(deconvlucy(lmobj.initimg, fspecial('gaussian',31,4.46))));
    imagesc(crop(img));
    axis equal; axis off;

end


function img2=crop(img)
img2 = img(end/4:end/4*3+1,end/4:end/4*3+1);
end

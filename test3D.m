% generate data
n=1e6;
data0 = rand(1,n) * 100 - 50; % x ~ -[50,50]
data0(2,:) = 0.1 * data0(1,:);
data0(3,:) = -log(rand(1,n));
data0(4,:) = -10;
I = data0(3,:) < 0.2;
data0(:,I) = [];
data0(3,:) = 5 ./ sqrt(data0(3,:));
data0(5,:) = data0(3,:) * 3;

tmp = data0; clear data0;

data0 = rand(1,n) * 100 - 50; 
data0(2,:) = -0.1 * data0(1,:);
data0(3,:) = -log(rand(1,n));
data0(4,:) = 10;
I = data0(3,:) < 0.2;
data0(:,I) = [];
data0(3,:) = 5 ./ sqrt(data0(3,:));
data0(5,:) = data0(3,:) * 3;

data0 = [data0 tmp];

data1 = data0;
data1(1:2,:) = data0(1:2,:) + randn(2, length(data0)).*data0(3,:);
data1(4,:) = data0(4,:) + randn(1,length(data0)).*data0(5,:);

lmobj=lmdatainit(data1(:,1:100:end),[2,4]);

% MAP estimator
img=lmdeconv(lmobj,40);
img0 = histimg(lmobj);

figure;
subplot(2,2,1);
imagesc(squeeze(sum(img0,2))); axis off equal;
subplot(2,2,2);
imagesc(squeeze(sum(img,2))); axis off equal;
subplot(2,2,3);
imagesc(squeeze(sum(img0,3))); axis off equal;
subplot(2,2,4);
imagesc(squeeze(sum(img,3))); axis off equal;
drawnow;

% MMSE estimator
disp ('Calculate MMSE estimator. This will take a while...')
s = lmsample(lmobj,1,[],1000);
st = s;
for i = 1:10
    sn = lmsample(lmobj,1000,s);
    st = st + sum(sn,4);
    s = sn(:,:,:,end);
    %disp('tick');
end

img = st;
figure;
subplot(2,2,1);
imagesc(squeeze(sum(img0,2))); axis off equal;
subplot(2,2,2);
imagesc(squeeze(sum(img,2))); axis off equal;
subplot(2,2,3);
imagesc(squeeze(sum(img0,3))); axis off equal;
subplot(2,2,4);
imagesc(squeeze(sum(img,3))); axis off equal;
drawnow;


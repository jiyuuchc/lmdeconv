function f=frc(img1,img2,nn)

if (~exist('nn', 'var'))
  nn=100;
end

fftimg1=fft2(img1)-mean(img1(:));
fftimg2=conj(fft2(img2)-mean(img2(:)));
fftc=fftimg1.*fftimg2;
ffta1=abs(fftimg1).^2;
ffta2=abs(fftimg2).^2;

[h,w]=size(fftc);
x=[0:(w-1)/2 fliplr(-1:-1:-(w/2))]/w;
y=[0:(h-1)/2 fliplr(-1:-1:-(h/2))]/h;
[X,Y]=meshgrid(x,y);
R=sqrt(X.^2+Y.^2);
[N,edges]=histcounts(R(:),nn);
Rd = discretize(R,edges);
for i = 1:nn
    frc(i)=sum(fftc(Rd==i))/sqrt(sum(ffta1(Rd==i)) * sum(ffta2(Rd==i)));
end
f(1,:)=(edges(1:end-1)+edges(2:end))/2;
f(2,:)=real(frc);

function [newParticles, fval]=particleIter(particles, pixelSize, nsamples, scale, prior)
% newParticles = particleIter(particles, pixelSize, nsamples, scale)
% Perform particle registration using SMLM data. This function can be
% called multiple times to iteratively improve the registration
% Input:
%   particles: cell array of particle data (Nx3 matrix);Discretized.
%   pixelSize : pixelsize used for theta samples. 
%   nsamples: number of samples drawn
%   scale: scale sigma by a factor, may improve convergence. default to 1.
% Output:
%   newParticles: New particle coordinates with improved registration

skipping = 5; % skip every 5 samples to improve mixing
burn_in = 2000; % take 2000 samples as burn in

if(~exist('prior','var'))
    prior = 1;
end

if (~exist('scale','var'))
    scale = 1.0;
end

alldata = cat(1, particles{:});
alldata(3,:)=alldata(3,:)*scale;
lmobj = lmdatainit(alldata', pixelSize);

disp('Sampling');
tic;
s = lmsample(lmobj,1,[],burn_in,prior);
l = zeros([size(s) length(lmobj.psfs)]);
for k = 1:1000:nsamples
    nn = min(1000,nsamples-k+1);
    samples = lmsample(lmobj, nn, s, skipping, prior);
    s = samples(:,:,end);
    
    tmp = zeros(size(samples));
    for j = 1:length(lmobj.psfs)
        psf = lmobj.psfs{j};
        for i = 1:size(samples,3)
            tmp(:,:,i) = filter2(psf, samples(:,:,i)+eps);
        end
        l(:,:,j) = l(:,:,j) + mean(log(tmp),3);
    end
end
toc;

disp('registering');
tic;
newParticles = {};
idx = 1;
for i = 1:length(particles)
    nlocs = length(particles{i});
    %HACKISH
    tmp = double(lmobj.data(:, idx:idx+nlocs-1));
    idx = idx + nlocs;
    [p, fval(i)] = fminunc(@(p) costfunc(p, tmp', l), [0 0 0]);
    t = affine2d([cos(p(3)) sin(p(3)) 0; -sin(p(3)) cos(p(3)) 0; p(1) p(2) 1]);
    tmp = double(particles{i});
    tmp(:,1:2)= t.transformPointsForward(tmp(:,1:2));
    newParticles{i} = tmp;
end

toc;

function ll = costfunc(p, data, l)
dt = p(3);
t = affine2d([cos(dt) sin(dt) 0; -sin(dt) cos(dt) 0; p(1) p(2) 1]);
d2 = t.transformPointsForward(data(:,1:2))+0.5;
ll = - sum(interp3(l, d2(:,1), d2(:,2), data(:,3), 'linear', log(eps)));



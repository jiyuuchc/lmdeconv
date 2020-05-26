function [newParticles, fval]=particleIter(particles, pixelSize, nsamples, target, scale, prior)
% newParticles = particleIter(particles, pixelSize, nsamples, target, scale)
% Perform particle registration using SMLM data. This function can be
% called multiple times to iteratively improve the registration
% Input:
%   particles: cell array of particle data (Nx3 matrix);Discretized.
%   pixelSize : pixelsize used for theta samples. 
%   nsamples: number of samples drawn
%   target: optional. if specified, register all particles against a single
%           target particle. Nx3 array.   
%   scale: scale sigma by a factor, may improve convergence. default to 1.
% Output:
%   newParticles: New particle coordinates with improved registration

skipping = 5; % skip every 5 samples to improve mixing
burn_in = 2000; % take 2000 samples as burn in

if(~exist('prior','var'))
    prior = .5;
end

if (~exist('scale','var'))
    scale = 1.0;
end

alldata = cat(1, particles{:});

if (exist('target','var') && ~isempty(target))
    alldata = [target; alldata];
end

alldata(3,:)=alldata(3,:)*scale;
lmobj = lmdatainit(alldata', pixelSize);
ds = lmobj.S(2)-lmobj.S(1);
sigmaedges = [lmobj.S - ds/2, inf];

if (exist('target','var') && ~isempty(target))
    lmobj.data  = lmobj.data(:,1:size(target,1)); %only keep target for sampling
end

disp('---------Expectation------');
tic;
disp('burn in...');
s = lmsample(lmobj,1,[],burn_in,prior);
disp('burn in...done');
disp('Sampling...');
l = zeros([size(s) length(lmobj.psfs)]);
for k = 1:1000:nsamples
    nn = min(1000,nsamples-k+1);
    samples = lmsample(lmobj, nn, s, skipping, prior);
    s = samples(:,:,end);
    
    tmp = zeros(size(samples));
    for j = 1:length(lmobj.psfs)
        psf = lmobj.psfs{j};
        for i = 1:size(samples,3)
            tmp(:,:,i) = filter2(psf, samples(:,:,i));
        end
        l(:,:,j) = l(:,:,j) + mean(log(tmp),3);
    end
end
disp('Sampling...done');
toc;

disp('---------Maximization----------');
disp('registering....');
tic;
newParticles = particles;
fval = zeros(1,length(particles));
options = optimoptions(@fminunc,'Display','none');
for i = 1:length(particles)
    curdata = particles{i};
    curdata(:,3) = discretize(curdata(:,3),sigmaedges);
    func = @(p) costfunc(p, curdata, l, lmobj.X, lmobj.Y, 1:size(l,3));
    [p, fval(i)] = fminunc(func, [0,0,0], options);
    t = affine2d([cos(p(3)) sin(p(3)) 0; -sin(p(3)) cos(p(3)) 0; p(1) p(2) 1]);
    newParticles{i}(:,1:2) = t.transformPointsForward(curdata(:,1:2));
    if (mod(i, 100) == 0)
        disp(['finished ', int2str(i), ' particles']);
    end
end
disp('registering...done');
toc;

function ll = costfunc(p, data, V, X, Y, Z)
dt = p(3);
t = affine2d([cos(dt) sin(dt) 0; -sin(dt) cos(dt) 0; p(1) p(2) 1]);
d2 = t.transformPointsForward(data(:,1:2));
ll = - sum(interp3(X, Y, Z, V, d2(:,1), d2(:,2), data(:,3), 'linear', log(eps(0))));

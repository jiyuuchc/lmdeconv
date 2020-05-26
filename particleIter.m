function [newParticles, fval]=particleIter(particles, pixelSize, nsamples, varargin)
% newParticles = particleIter(particles, pixelSize, nsamples, ...)
% Perform particle registration using SMLM data. This function can be
% called multiple times to iteratively improve the registration
% Input:
%   particles: cell array of particle data (Nx3 matrix);Discretized.
%   pixelSize : pixelsize used for theta samples. 
%   nsamples: number of samples drawn
%   Optional Parameters:
%       'Target': If specified, register all particles against a single
%           target particle. Nx3 array.
%       'Scale': scale sigma by a factor, may improve convergence. default to 1.
%       'BurnIn': # of samples as burn-in, default 2000
%       'Skipping': # of skipped samples (to decrease correlation). default 5
% Output:
%   newParticles: New particle coordinates with improved registration

parser = inputParser;
isnumericposscalar = @(p) isnumeric(p) && isscalar(p) && p > 0;
addRequired(parser, 'particles', @(p) iscell(p));
addRequired(parser, 'pixelSize', isnumericposscalar);
addRequired(parser, 'nsamples', isnumericposscalar);
addParameter(parser, 'Scale', 1, isnumericposscalar);
addParameter(parser, 'Prior', 0.5, isnumericposscalar);
addParameter(parser, 'BurnIn', 2000,isnumericposscalar);
addParameter(parser, 'Skipping', 5, isnumericposscalar);
addParameter(parser, 'Target', [], @(p) isnumeric(p) && size(p,2) == 3);
addParameter(parser, 'FitFun', @fit_gradient);
parse(parser, particles, pixelSize, nsamples, varargin{:});

skipping = parser.Results.Skipping;
burn_in = parser.Results.BurnIn;
prior = parser.Results.Prior;
scale = parser.Results.Scale;
target = parser.Results.Target;

alldata = cat(1, particles{:});

if (~isempty(target))
    alldata = [target; alldata];
end

alldata(3,:)=alldata(3,:)*scale;
lmobj = lmdatainit(alldata', pixelSize);

if (~isempty(target))
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
for i = 1:length(particles)
    curdata = particles{i};
    [t, fval(i)] = parser.Results.FitFun(curdata, l, lmobj);
    newParticles{i}(:,1:2) = t.transformPointsForward(curdata(:,1:2));
    if (mod(i, 100) == 0)
        disp(['finished ', int2str(i), ' particles']);
    end
end
disp('registering...done');
toc;

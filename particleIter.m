function [newParticles, model, fval]=particleIter(particles, pixelSize, nsamples, varargin)
% newParticles = particleIter(particles, pixelSize, nsamples, ...)
% newParticles, model = particleIter(particles, pixelSize, nsamples, model, ...)
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
addParameter(parser, "Model", 1, @(p) isvector(p) || isnumericposscalar(p));
addParameter(parser, 'Scale', 1, isnumericposscalar);
addParameter(parser, 'Prior', 0.5, isnumericposscalar);
addParameter(parser, 'BurnIn', 2000,isnumericposscalar);
addParameter(parser, 'Skipping', 5, isnumericposscalar);
addParameter(parser, 'Target', [], @(p) ismatrix(p) && size(p,2) == 3);
addParameter(parser, 'FitFun', @fit_gradient);
parse(parser, particles, pixelSize, nsamples, varargin{:});

skipping = parser.Results.Skipping;
burn_in = parser.Results.BurnIn;
prior = parser.Results.Prior;
scale = parser.Results.Scale;
target = parser.Results.Target;

if (isscalar(parser.Results.Model))
    nModels = max(1, floor(parser.Results.Model));
    model = randi(nModels, 1, length(particles));
else
    model = floor(parser.Results.Model);
    if length(model) ~= length(particles)
        error('length of model assignement should match particle data');
    end
    if (min(model) < 1)
        error('model assignment should be >= 1');
    end
    nModels = max(model);
end

disp('---------Expectation------');
alldata = cat(1, particles{:});
lmobj_all = lmdatainit(alldata, pixelSize);
template = cell(1, nModels);
for m = 1:nModels
    collections = model == m;
    alldata = cat(1, particles{collections});

    if (~isempty(target) && nModels == 1)
        alldata = [target; alldata];
    end

    alldata(3,:)=alldata(3,:)*scale;
    obj = lmdatainit(alldata, lmobj_all);

    if (~isempty(target) && nModels == 1)
        obj.data  = obj.data(:,1:size(target,1)); %only keep target for sampling
    end
    
    if (nModels > 1)
        disp(['Model #' int2str(m)]);
    end
    
    tic;
    disp('burn in...');
    s = lmsample(obj,1,'skip',burn_in,'Prior',prior);
    disp('burn in...done');
    disp('Sampling...');
    l = zeros([size(s) length(obj.psfs)]);
    for k = 1:1000:nsamples
        nn = min(1000,nsamples-k+1);
        samples = lmsample(obj, nn, s, 'Skip',skipping, 'Prior',prior, 'Quite', true);
        s = samples(:,:,end);

        tmp = zeros(size(samples));
        for j = 1:length(obj.psfs)
            psf = obj.psfs{j};
            for i = 1:size(samples,3)
                tmp(:,:,i) = filter2(psf, samples(:,:,i));
            end
            l(:,:,j) = l(:,:,j) + mean(log(tmp),3);
        end
        disp(['Drawn ' int2str(k+999) ' samples.']);
    end
    template{m} = l;
    disp('Sampling...done');
    toc;
end

disp('---------Maximization----------');
disp('registering....');
tic;
newParticles = particles;
fval = zeros(1,length(particles));
for i = 1:length(particles)
    curdata = particles{i};
    vmin = realmax;
    for m = 1: nModels
        [t, v] = parser.Results.FitFun(curdata, template{m}, lmobj_all);
        if (v < vmin)
            vmin = v;
            tbest = t;
            mbest = m;
        end
    end
    fval(i) = vmin; 
    newParticles{i}(:,1:2) = tbest.transformPointsForward(curdata(:,1:2));
    model(i) = mbest;
    if (mod(i, 100) == 0)
        disp(['finished ', int2str(i), ' particles']);
    end
end
disp('registering...done');
toc;

function [newParticles, model, fval]=particleIter(particles, pixelSize, nsamples, varargin)
% newParticles = particleIter(particles, pixelSize, nsamples, ...)
% [newParticles, model] = particleIter(particles, pixelSize, nsamples, 'Model', model, ...)
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
%       'FitFun': optimization routine: @fit_gradient (default), @fit_ps, @fit_global
%       'Model': either a scalar (number of models) or a vector (model assignments)
% Output:
%   newParticles: New particle coordinates with improved registration
%   model: model assignments

parser = inputParser;
parser.KeepUnmatched = true;
isnumericposscalar = @(p) isnumeric(p) && isscalar(p) && p > 0;
addRequired(parser, 'particles', @(p) iscell(p));
addRequired(parser, 'pixelSize', isnumericposscalar);
addRequired(parser, 'nsamples', isnumericposscalar);
addParameter(parser, "Model", 1, @(p) isvector(p) || isnumericposscalar(p));
addParameter(parser, 'Scale', 1, isnumericposscalar);

addParameter(parser, 'Target', [], @(p) ismatrix(p) && size(p,2) == 3);
addParameter(parser, 'FitFun', @fit_gradient);
parse(parser, particles, pixelSize, nsamples, varargin{:});

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

    alldata(3,:)=alldata(3,:) * parser.Results.Scale;
    obj = lmdatainit(alldata, lmobj_all);

    if (~isempty(target) && nModels == 1)
        obj.data  = obj.data(:,1:size(target,1)); %only keep target for sampling
    end

    if (nModels > 1)
        disp(['Model #' int2str(m)]);
    end

    tic;
    template{m} = intloglikelihood(obj, nsamples, parser.Unmatched);
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

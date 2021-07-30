function img = mmseimg(lmobj, varargin)
% img = mmseimg(lmobj, nsamples, ...)
% Perform statitically accurate rendering of SMLM data using MMSE estimator
%
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   nsamples - # of iterations, optional, default 20000
%   Optional parameters:
%     'BurnIn' - # of samples for burn in
%     'BatchSize' - Batchsize for sampling. Limited by memory.
%     'Skipping' - # of samples drawn for every sample kept. Default 1(not skipping)
%     'Prior' - Dirichlet prior, Optional.
% Outputs -
%   img - MMSE rendering

isPosInt = @(p) validateattributes(p, {'numeric'},{'scalar','positive','integer'});
isnumericposscalar = @(p) isnumeric(p) && isscalar(p) && p > 0;

parser = inputParser;
addRequired(parser, 'lmobj', @(p) isstruct(p));
addOptional(parser, 'nsamples', 20000, isPosInt);

addParameter(parser, 'Skipping', 1, isnumericposscalar);
addParameter(parser, 'BurnIn', 5000, isnumericposscalar);
addParameter(parser, 'Prior', 0.5, @(p) isnumeric(p) && (isscalar(p) || ismatrix(p)));
addParameter(parser, 'BatchSize', 1000, isPosInt);
addParameter(parser, 'Quiet', false, @(p) islogical(p));

skipping = parser.Results.Skipping;
burn_in = parser.Results.BurnIn;
prior = parser.Results.Prior;
batchSize = parser.Results.BatchSize;

disp('burn in...');
s = lmsample(lmobj,1,'skip',burn_in,'Prior',prior);
img = zeros(size(s));

disp('Sampling...');
for k = 1:batchSize:nsamples
    nn = min(batchSize, nsamples-k+1);
    samples = lmsample(lmobj, nn, s, 'Skip',skipping, 'Prior',prior, 'Quiet', true);
    s = samples(:,:,end);

    img = img + sum(samples, 3);
    if (~ parser.Results.Quiet)
      disp(['Drawn ' int2str(k + batchSize - 1) ' samples.']);
    end
end
img = img / nsamples;

function model = intloglikelihood(lmobj, nsamples, varargin)
% compute the integrated (expectation) loglikelihood under various gaussian filtering

isPosInt = @(p) validateattributes(p, {'numeric'},{'scalar','positive','integer'});
isnumericposscalar = @(p) isnumeric(p) && isscalar(p) && p > 0;

parser = inputParser;
addRequired(parser, 'lmobj', @(p) isstruct(p));
addRequired(parser, 'nsamples', isnumericposscalar);
addParameter(parser, 'Skipping', 5, isnumericposscalar);
addParameter(parser, 'BurnIn', 2000, isnumericposscalar);
addParameter(parser, 'Prior', 0.5, @(p) isnumeric(p) && (isscalar(p) || ismatrix(p)));
addParameter(parser, 'BatchSize', 500, isPosInt);
addParameter(parser, 'UseGpu', false, @(p) islogical(p));
addParameter(parser, 'Quiet', false, @(p) islogical(p));

parse(parser, lmobj, nsamples, varargin{:});

skipping = parser.Results.Skipping;
burn_in = parser.Results.BurnIn;
prior = parser.Results.Prior;
batchSize = parser.Results.BatchSize;
useGpu = parser.Results.UseGpu;

disp('burn in...');
s = lmsample(lmobj,1,'skip',burn_in,'Prior',prior);

disp('Sampling...');
if (useGpu)
  model = zeros([size(s) length(lmobj.S)], 'gpuArray');
else
  model = zeros([size(s) length(lmobj.S)]);
end
for k = 1:batchSize:nsamples
    nn = min(batchSize, nsamples-k+1);
    samples = lmsample(lmobj, nn, s, 'Skip',skipping, 'Prior',prior, 'Quiet', true);
    s = samples(:,:,end);

    if (useGpu)
      samples = gpuArray(samples);
    end
    for j = 1:length(lmobj.S)
        psfsize = lmobj.S(j) / lmobj.pixelsize;
        model(:,:,j) = model(:,:,j) + sum(log(imgaussfilt(samples, psfsize)), 3);
    end

    if (~ parser.Results.Quiet)
      disp(['Drawn ' int2str(k + batchSize - 1) ' samples.']);
    end
end

%    if (useGpu)
%      model = gather(model);
%    end
model = model / nsamples;

function img = lmdeconv(lmobj, varargin)
% img = lmdeconv(locdata, iters, prev, prior)
% Perform statitically accurate rendering of SMLM data using ML estimator
% 
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iter - # of iterations
%   prev - result from previos call of this function. Optional
%   prior - Dirichlet prior, Optional.
% Outputs -
%   img - deconv img

parser = inputParser;
isPosInt = @(p) validateattributes(p, {'numeric'},{'scalar','positive','integer'});
addRequired(parser, 'lmobj', @(p) isstruct(p));
addOptional(parser, 'iters', 100, isPosInt);
addParameter(parser, 'Init', [], @(p) validateattributes(p, {'numeric'},{'2d'}));
addParameter(parser, 'Prior', .5, @(p) isnumeric(p) && (isscalar(p) || ismatrix(p)));
addParameter(parser, 'Quiet', false, @(p) islogical(p));
parse(parser, lmobj, varargin{:});

iters = parser.Results.iters;
prev = parser.Results.Init;
prior = parser.Results.Prior;
if(isempty(prev))
    prev = ones(lmobj.imgsize);
end

for i = 1:iters
    if (length(lmobj.imgsize) == 2)
        prev = lmdeconvmex(lmobj.data, prev, lmobj.psfs, double(prior));
    else
        prev = lmdeconv3dmex(lmobj.data, prev, lmobj.psfs, lmobj.zpsfs, double(prior));
    end
    if (mod(i,50) == 0 && ~ parser.Results.Quiet)
        disp(['finished ' int2str(i) ' iterations']);
    end
end

img = prev;

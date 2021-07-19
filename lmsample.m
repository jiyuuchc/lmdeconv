function samples = lmsample(lmobj, iters, varargin)
% samples = lmdsample(locdata, iters, prev, ...)
% Sample the posterior distribution on theta based on SMLM data
%
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iter - # of samples to be drawn
%   prev - Optional. result from previos iteration. 
%   Optional parameters:
%       'Skip' - Optional. skipping samples, default 1
%       'Prior' - Dir prior. can be a scaler or a matrix.
%
% Outputs -
%   samples - samples

parser = inputParser;
isPosInt = @(p) validateattributes(p, {'numeric'},{'scalar','positive','integer'});
addRequired(parser, 'lmobj', @(p) isstruct(p));
addRequired(parser, 'iters', isPosInt);
addOptional(parser, 'prev', [], @(p) validateattributes(p, {'numeric'},{'2d'}));
addParameter(parser, 'Skip', 1, isPosInt);
addParameter(parser, 'Prior', .5, @(p) isnumeric(p) && (isscalar(p) || ismatrix(p)));
addParameter(parser, 'Quiet', false, @(p) islogical(p));
parse(parser, lmobj, iters, varargin{:});

prev = parser.Results.prev;
prior = parser.Results.Prior;
skip = parser.Results.Skip;

if(isempty(prev))
    prev = histimg(lmobj);
else
    if ~ all(size(prev) == lmobj.imgsize)
        error('sample size does not match lmobj data');
    end
end

samplesize = [size(prev) iters];
samples = zeros(samplesize);
for i = 1:iters * skip
    if (length(lmobj.imgsize) == 2)
        prev = lmsamplermex(lmobj.data, prev, lmobj.psfs, double(prior));
        if (mod(i,skip)==0)
            samples(:,:,i/skip) = prev;
        end
    else
        prev = lmsampler3dmex(lmobj.data, prev, lmobj.psfs, lmobj.zpsfs, double(prior));
        if (mod(i,skip)==0)
            samples(:,:,:,i/skip) = prev;
        end
    end
    if (mod(i,1000 * skip) == 0 && ~parser.Results.Quiet)
        disp(['finished ' int2str(i/skip) ' iterations']);
    end
end

function img = mapimg(lmobj, varargin)
% img = mapimg(lmobj, iters, ...)
% Perform statitically accurate rendering of SMLM data using MAP estimator
%
% Inputs -
%   lmobj - Object returned by the lmdatainit() function
%   iters - # of iterations, optional, default 100
%   Optional parameters:
%     'Init' - result from previos call of this function. Optional
%     'Prior' - Dirichlet prior, Optional.
% Outputs -
%   img - MAP rendering

img = lmdeconv(lmobj, varargin{:});

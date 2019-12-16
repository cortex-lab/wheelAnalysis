function ret = gausswin(Nech, a)
% WHEEL.GAUSSWIN(NECH, A) Gaussian window
%   Equivalent behaviour as the sigproc function.  Returns an N-point
%   Gaussian window with a given alpha value.
%
%   Inputs:
%     Nech - The number of points in the window
%     a - The alpha value, where alpha is the reciprocal of the standard
%     deviation; a measure of the width of the window's Fourier transform.
%     Default = 2.5.  The reulting s.d. = (Nech-1) / (2*a)
%
%   Example:
%     win = wheel.gausswin(5);
%     >> [0.0439 0.4578 1.0000 0.4578 0.0439]
if nargin <=1, a=2.5; end
x = linspace(-1,1,Nech)';
ret = exp(-0.5 .* (a.* x).^2);
function ret = gausswin(Nech, a)
% wheelmoves.gausswin()
% equivalent behaviour as the sigproc function
if nargin <=1, a=2.5; end
x = linspace(-1,1,Nech)';
ret = exp(-0.5 .* (a.* x).^2) ;
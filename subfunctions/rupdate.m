function [r F] = rupdate(varargin)
% [r F]= rupdate(x,n,v,tau,'L1')
% [r F]= rupdate(x,n,v,tau,'L2')

x = varargin{1};
n = varargin{2};
v = varargin{3};
tau = varargin{4};
if strcmp(varargin{5},'L1')
   abar = sum(acos(v'*x/norm(v)))/(n);
    if abs(pi/2 - abar) - tau > 0
        r = pi/2 - sign(pi/2-abar) * (abs(pi/2 - abar) - tau);
    else
        r = pi/2;
    end
    F = mean((acos(v'*x/norm(v)) - r).^2)/2 + tau  * abs(pi/2-r) ;
else  % then L2
    r = (sum(acos(v'*x/norm(v)))/(n)  + pi/2*tau)/(1+tau);
    F = mean((acos(v'*x/norm(v)) - r).^2)/2 + tau/2 * (pi/2-r)^2;
end
function rr = LMFsphereResPen(center)
% LMFSPHERERES : Auxiliary function for LMFsphereFit.
%              -Calculates residuals of circle fit for the given center.
%
%   See also LMFsphereFit, LMFnlsq.
% Last updated Aug 10, 2012
% Sungkyu Jung

global tpdata
global Radius
n = size(tpdata,2);
di = sqrt(sum((tpdata - repmat(center,1,n)).^2));
rr = (di-Radius)';


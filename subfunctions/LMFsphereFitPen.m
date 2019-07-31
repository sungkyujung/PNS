function [center] = LMFsphereFitPen(x,c0,r,loss,TOL)
global tpdata Radius;

[d n] = size(x);
c0 = c0/norm(c0); % normalize the new candidate
rot = rotMat(c0);  % rotation matrix : c0 -> North Pole
if strcmp(loss,'n');
    tpdata = NaiveExtrinsicLog(rot*x);
    Radius = 2*sin(r/2) ;
else
    tpdata = LogNPd(rot*x); % Tangent projection by Log map
    Radius = r;
end

newCenterTp=LMFnlsq('LMFsphereResPen',zeros(d-1,1),'Display',0,'MaxIter',50,'XTol',TOL);

if strcmp(loss,'n');
    newCenter = NaiveExtrinsicExp(newCenterTp);
else
    newCenter = ExpNPd(newCenterTp); % Bring back to the sphere by Exp map
end

center = rot\newCenter;  % (inv(rot)*newCenter) rotate back the newCenter

clear data;
clear Radius;

function [center, r] = getSubSpherePen(data,tau,loss,penalty,maxcnt,TOL)

if nargin == 4
    maxcnt = 100;
    TOL = 1e-9;
end

d= size(data,1);

% use Slicing Loss to get initial values
[c01] = getSubSphereSlicing(data,tau,'null',10,1e-5);
[c02] = getSubSphereSlicing(data,tau,'G',10,1e-5);
[c03] = getSubSphereSlicing(data,tau,penalty,10,1e-5);
c0vec = [c01 c02 c03];

cmpinitials = zeros(d+2,3);
for init = 1:3
    c0 = c0vec(:,init);
    cnt =0;
    [r, Fnow]= rupdate(data,c0,tau,loss,penalty);
    err = 1;
    %rvec = []; 
    while err > TOL
        % Given the r, update center;
        newcenter = LMFsphereFitPen(data,c0,r,loss,TOL) ;
        % Given the center, update r
        [r, Fnext]= rupdate(data,newcenter,tau,loss,penalty);
        err = abs(Fnow - Fnext) ;
        Fnow = Fnext;
        c0 = newcenter;
        cnt = cnt+1;
        %rvec = [rvec r];
        if cnt > maxcnt;
            disp(['Message from getSubSpherePen.m:   iteration reached its maximum size with error ' num2str(err)]);
            disp(['init = ' num2str(init)])
            break;
        end
    end
    cmpinitials(1:d,init) = newcenter;
    cmpinitials(d+1,init) = r;
    cmpinitials(d+2,init) = Fnow;
end

[tmp mininit]=min(cmpinitials(d+2,:));
center = cmpinitials(1:d,mininit);
r = cmpinitials(d+1,mininit);




if r > pi/2;
    center = -center;
    r = pi - r;
end % adjust radius






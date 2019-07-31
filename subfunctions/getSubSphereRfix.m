function [center, F] = getSubSphereRfix(data,r,loss)
%   loss
%            'i'  (default) Intrinsic squared loss
%            'n'  Naive extrinsic squared loss
%            's'  Slicing squared loss
if nargin==0 && nargout==0, help getSubSphereRfix, return, end     %   Display help

[d n] = size(data);
TOL = 1e-9;
maxcnt = 100;




c0vec = zeros(d,2);

% The first initial value is from 'null'
cdata = data - repmat(mean(data,2),1,n);
[U dd] = svd(cdata);
[tmp minindex]=min(diag(dd));
c0vec(:,1) = U(:,minindex);

% The second initial value is from 'G'
[U dd] = svd(data);
[tmp minindex]=min(diag(dd));
c0vec(:,2) = U(:,minindex);

switch loss
    case 's'
        b = cos(r);
        S = data*data' / n;
        xbar = mean(data,2);
        [U L]=eig(S);
        lambdaS = min(diag(L));
        xstar = U'*xbar;
        lend = -b^2*sum(xbar.^2);
        rend = lambdaS - TOL;
        if abs(b) < TOL;
            center =c0vec(:,2);
        else
            diff = 1;
            while diff > TOL;
                mid = (lend+rend) / 2;
                fmid = b^2* sum((xstar.^2) ./ (diag(L) - mid).^2) - 1;
                if fmid <= 0
                    diff = mid - lend;
                    lend = mid;
                else
                    diff = rend - mid;
                    rend = mid;
                end
            end
            center =   b*((S - mid*eye(d))\xbar);
        end
    case 'i'
        
        cmpinitials = zeros(d+2,2);
        for init = 1:2
            c0 = c0vec(:,init);
            b = r;
             
            cnt =0;
            err = 1;
            
            Fcurrent = mean(   (acos(c0'*data/norm(c0))        - b).^2  )/2;
            
            while err > TOL
                % Given r, update center
                newcenter = LMFsphereFitPen(data,c0,r,loss,TOL) ;
                Fnew = mean(   (acos(newcenter'*data/norm(newcenter))        - b).^2  )/2;
                err = abs(Fcurrent - Fnew) ;
                Fcurrent = Fnew;
                c0 = newcenter;
                cnt = cnt+1;
                if cnt > maxcnt;
                    break;
                end
            end
            cmpinitials(1:d,init) = newcenter;
            cmpinitials(d+1,init) = r;
            cmpinitials(d+2,init) = Fnew;
        end
        
        [tmp mininit]=min(cmpinitials(d+2,:));
        center = cmpinitials(1:d,mininit);
        
    case 'n'
        cmpinitials = zeros(d+2,3);
        for init = 1:2
            c0 = c0vec(:,init);
            b = r;
             
            cnt =0;
            err = 1;
            Fcurrent = mean(  (sqrt( 2 - 2*(c0'*data/norm(c0))  ) - b).^2  )/2;
            
            while err > TOL
                % Given r, update center
                newcenter = LMFsphereFitPen(data,c0,r,loss,TOL) ;
                Fnew = mean(  (sqrt( 2 - 2*(newcenter'*data/norm(newcenter))  ) - b).^2  )/2;
                err = abs(Fcurrent - Fnew) ;
                Fcurrent = Fnew;
                c0 = newcenter;
                cnt = cnt+1;
                if cnt > maxcnt;
                    break;
                end
            end
            cmpinitials(1:d,init) = newcenter;
            cmpinitials(d+1,init) = r;
            cmpinitials(d+2,init) = Fnew;
        end
        
        [tmp mininit]=min(cmpinitials(d+2,:));
        center = cmpinitials(1:d,mininit);
end



F = Fcurrent;





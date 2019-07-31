function [Tarray, distancefromtheaxis,NSOrthaxis,NSradius]= PNSanalysisofFit(resmat,PNS,k)
% PNSanalysisofFit 
% [Tarray, distancefromtheaxis,NSOrthaxis] = PNSe2s(data,PNS)
%      where 'data' is d x m data matrix in PNS coordinate system (for any
%      m >= 1), 'PNS' is the structural array from PNSmain.m
%       k < d represents the number of consecutive dimension where you want 
%              to investigate (default k = d - 1). For example, 
%       k = 2 (for S^3, i.e., d = 3), you get projected data on S2 and S3. 
%       k = 4 (for S^9, i.e., d = 9), you get projected data on S2, S3, S4 and S5.
 
[dm, n] = size(resmat);
% dm is the intrinsic dimension of the sphere
%    or the reduced sphere in HDLSS case.
% n  is the sample size.

if nargin < 3; 
    k = dm-1;
end
if k < 1; k = 1; end

NSOrthaxis = flipud(PNS.orthaxis(1:end-1));
geodmean = PNS.orthaxis{end};

% "standardize" the coordinates
res = resmat./repmat(flipud(PNS.radii),1,n);
% res(k,i) is the deviance of ith sample from fitting (k-1)-subsphere
% from S^k. For example, res(2,1) is the (signed) deviance of the
% (i)th observation by fitting an actual, possibly small, circle from
% projected data on S^2. We do not utilize the first row of 'res'.
if dm > 2;
    NSradius = flipud(PNS.dist);
else
    NSradius = PNS.dist;
end

distancefromtheaxis = repmat(NSradius,1,n) + res(2:end,:);

% figure(2);clf;
% histJ(distancefromtheaxis(1,:))

% iteratively mapping back to S^d

% to S^1
if dm > 0; 
    T = [cos(geodmean+res(1,:)) ; sin(geodmean+res(1,:))];
end

% S^1 to S^2
if dm > 1
T = rotMat(NSOrthaxis{1})'* ...
    [repmat(sin(NSradius(1)+res(2,:)),2,1).* T ;
    cos(NSradius(1)+res(2,:))];
end
Tarray = cell(k,1);
Tarray{1} = T;


% S^2 to S^d
if dm > 2
for i=1:dm-2;
    T = rotMat(NSOrthaxis{i+1})'*...
        [repmat(sin(NSradius(i+1)+res(i+2,:)),2+i,1).*T ;
        cos(NSradius(i+1)+res(i+2,:))];
    if i < k
    Tarray{i+1} = T;
    end
end
end






% 
% if isempty(PNS.basisu) == 0
%     % Then this is the HDLSS case
%     T = PNS.basisu*T;
% end 
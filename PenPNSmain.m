function [resmat PNS]= PenPNSmain(varargin)
% PENPNSMAIN Analysis of Penalized Principal Nested Spheres for directional
% data with penalty term 
% [resmat PNS]= PenPNSmain(data)    
% [resmat PNS]= PenPNSmain(data,tau) 
% [resmat PNS]= PenPNSmain(data,tau,'option') 
%
% Input:
%   data     (d+1) x n data matrix where each column is a unit vector.
% 
%   tau     1 (default) : tuning parameter
% 
%   option   'L2' L2 penalty 
%            'L1' L1 penalty (default)
%
% Output:
%   resmat  The commensurate residual matrix (X_PNS). Each entry in row k
%           works like the kth principal component score.
%
%   PNS     Matlab structure array with following fields:
%       mean     : location of the PNSmean
%       radii    : size (radius) of PNS
%       orthoaxis: orthogonal axis 'v_i' of subspheres
%       dist     : distance 'r_i' of subspheres
%       pvalues  : p-values of LRT and parametric boostrap tests (if any)
%       ratio    : estimated ratios
%       itype    : type of methods for fitting subspheres (if applicable)
data = varargin{1};

if nargin == 1;
    tau = 0;   
    option = 'L1'; % L1
elseif nargin < 3
    tau = varargin{2};
    option = 'L1';  % L1
else 
    tau = varargin{2};
    option = varargin{3};     
end

[k n] = size(data); % data on (k-1)-sphere in Real^k
[uu ll]=svd(data);
maxd = find(diag(ll) < 1e-15,1);
if isempty(maxd) || k > n
    maxd = min(k,n)+1;
end
nullspdim = k - maxd + 1; % dimension of subspace that contains no data

d = k - 1; % dimension of the sphere
disp(['Message from PNSmain.m; dataset is on ' num2str(d) '-sphere.']);

if nullspdim > 0
    disp([' .. found null space of dimension ' num2str(nullspdim)...
        ',to be trivially reduced.']);
end
resmat = zeros(d,n); % d dimensional residual matrix
% there will be d-1 subspheres
orthaxis = cell(d-1,1);
dist = zeros(d-1,1);
pvalues = zeros(d-1,2);
ratio = zeros(d-1,1);

% (HDLSS case) fit nested great spheres for dimension reduction
% where no residual is present.
currentSphere = data;
for i = 1:nullspdim;
    oaxis = uu(:,end-i+1); % orthogonal axis for subsphere
    r = pi/2;              % distance for subsphere
    pvalues(i,:) = [NaN,NaN];        % No test performed (p-value is Not-a-Number)
    res = acos(oaxis'*currentSphere)-r; % residuals
    % save subsphere parameters
    orthaxis{i} = oaxis; dist(i) = r;
    % save residuals
    resmat(i,:) = res;
    % projection to subsphere and transformation to isomorphic sphere
    NestedSphere = rotMat(oaxis)*currentSphere;
    currentSphere = NestedSphere(1:(k-i),:)./...
        repmat(sqrt(1-NestedSphere(end,:).^2),k-i,1);
    % transform singular vectors accordingly
    uu = rotMat(oaxis)*uu;
    uu = uu(1:(k-i),:)./repmat(sqrt(1-uu(end,:).^2),k-i,1); 
end 

% Now do sequential dimension reduction
for i=nullspdim+1:(d-1)
    if i == 3
        
    end
    % estimate the best fitting subsphere with tuning parameter tau
    [center, r] = getSubSpherePen(currentSphere,tau,option) ; 
    
    if  r> pi/2 ; 
        r = pi - r;
        center = -center;
    end
    
    res = acos(center'*currentSphere)-r;
    % save subsphere parameters
    orthaxis{i} = center; dist(i) = r;
    % save residuals
    resmat(i,:) = res;
    % projection to subsphere and transformation to isomorphic
    % sphere
    NestedSphere = rotMat(center)*currentSphere;
    currentSphere = NestedSphere(1:(k-i),:)./...
        repmat(sqrt(1-NestedSphere(end,:).^2),k-i,1);
end

% currentSphere has (intrinsic) dimension 1
% compute PNSmean and deviations.

% parametrize 1-sphere to angles
S1toRadian = atan2(currentSphere(2,:),currentSphere(1,:));
% Geodesic mean of angles
meantheta = geodmeanS1(S1toRadian');
orthaxis{d} = meantheta;
% save deviations from PNSmean
resmat(d,:) =mod(S1toRadian - meantheta + pi,2*pi)-pi;


radii =1;
for i = 1:(d-1)
    radii = [radii; prod(sin(dist(1:i)))];
end
resmat = flipud(repmat(radii,1,n).*resmat); % scaled residuals

PNS.radii = radii;         % size (radius) of nested spheres from largest to smallest

PNS.orthaxis = orthaxis;    % orthogonal axis of (d-1) subspheres and the anglemean for PNSmean
PNS.dist = dist;            % distances for (d-1) subspheres
%PNS.pvalues = pvalues;      % d-1 pvalues from sequential tests
PNS.ratio = ratio;   % d-1 ratios estimated
PNS.basisu = [];
PNS.mean = PNSe2s(zeros(d,1),PNS); % PNSmean of the data

end

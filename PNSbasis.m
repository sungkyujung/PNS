function Vstar = PNSbasis(PNS)
% PNSBASIS : Obtain basis from the result of PNS.
% 
% V = PNSbasis(PNS) for PNS, the output from PNSmain.m, returns d x
%   min(d,n) matrix consisting of basis column vectors.
%
% This code is based on Jung, Dryden and Marron (2012), 'Analysis of
% Principal Nested Spheres' Biometrika, 99(3), 551-568. In particular,
% Lemma 3 in page 13 of Supplementary material provides the formula for 
% the computation of the basis.
%
% Sungkyu Jung, Feb 2013.
  
 
r = PNS.dist;          % collection of r_1,...,r_{d-1}
sinr = [sin(r) ; 1];
d = length(sinr);
V = PNS.orthaxis;      % collection of v_1,...,v_{d}
V{d+1} = [-sin(PNS.orthaxis{end}) ; cos(PNS.orthaxis{end})];  % first component;
V{d} = [cos(PNS.orthaxis{end}) ; sin(PNS.orthaxis{end})];     % mean 

% First obtain all v^dagger
d = length(sinr);
Vdagger = [V{d+1} V{d}]; % matrix of v^dagger that will grow inside the loop from 2 to d+1;
for j = d:-1:2;
    Vdagger = rotMat(V{j-1})'*[sinr(j-1) * Vdagger ; repmat(sqrt(1-(sinr(j-1)).^2),1,d-j+2)];
    Vdagger = [Vdagger V{j-1}];
end
Vdagger = fliplr(Vdagger);
% switch order, so that v1^dagger is at Vdagger(:,1);
% sum(Vdagger.^2)

cosr = [(cos(r)); 0];
% cosr.^2 + sinr.^2
prodsinr = [PNS.radii(1:end) ; PNS.radii(end)];

% Now compute V^star
Vstar = zeros(d+1);
scaledVstar = zeros(d+1);
for j = 1:(d+1)
    vv = Vdagger(:,j); 
    scaledVstar(:,j) = vv - scaledVstar(:,1:j-1)*cosr(1:j-1); 
    Vstar(:,j) = scaledVstar(:,j)/prodsinr(j);
end
% Vstar'*Vstar

Vstar = fliplr(Vstar); 
% switch order again so that leading columns contain the `most important'
% basis vectors.
 
if isempty(PNS.basisu) == 0
    % Then this is the HDLSS case
    Vstar = PNS.basisu*Vstar;
end
 

% % simulated example;
% n =50;
% theta = linspace(0,pi,n);
% data = 5*[cos(theta); sin(theta)] + randn(2,n); data = data/10;
% data = rotMat([1 2 0]',[0 0 1]')\ExpNPd(data);
% 
% [resmat PNS]=PNSmain(data);
% V = PNSbasis(PNS);
% projections = data'*V;
% 
% figure(1);clf;
% PNSgraphics(resmat,PNS);
% figure(2);clf;
% subplot(2,2,1);
% scatterS2(data);
% hold on;
% plot3([0 V(1,1)],[0 V(2,1)],[0 V(3,1)],'b','Linewidth',3)  
% plot3([0 V(1,2)],[0 V(2,2)],[0 V(3,2)],'b','Linewidth',2)   
% plot3([0 V(1,3)],[0 V(2,3)],[0 V(3,3)],'b','Linewidth',1)   
% title('data + basis vectors');
% subplot(2,2,2);
% scatter(projections(:,1),projections(:,2));axis equal;
% xlabel('on v1')
% ylabel('on v2')
% subplot(2,2,4);
% scatter(projections(:,1),projections(:,3));axis equal;
% xlabel('on v1')
% ylabel('on v3')
% subplot(2,2,3);
% scatter(projections(:,2),projections(:,3));axis equal;
% xlabel('on v2')
% ylabel('on v3')
function pvalue = vMFtest(data,R)
% VMFTEST tests if the underlying distribution is an isotropic distribution
% with a single mode, assuming von-Mises Fisher distribution. 
% pvalue = vMFtest(data)
% pvalue = vMFtest(data,R), with d x n data matrix, give pvalue from a
% parametric boostrap with R (default 100) sets of random samples generated
% from the null distribution;
%   H0 : data ~ von Mises Fisher distribution (mu, kappa)
%   H1 : not H0
%
% using randvonMisesFisherm.m, getSubsphere.m
%
%-------------------------------------------------------
% NOTE:
% If you have Parallel Computing Toolbox, you can use "parfor" 
% by typing "matlabpool open" in Command Window before runnning scripts
% such as PNSmain.m and PNSshape.m.
% When finished learning PNS, type "matlabpool close" in Command Window.
% 
% I tried to insert "matlabpool open" in script files, but it didn't work.
%
% Takayuki Iwamoto(tiwamoto@cs.cmu.edu)
% 11/6/2013


[d n] = size(data);
if nargin < 2;
    R = 100;
end

global g_use_parallel;

% 1. test statistics
sumx = sum(data,2); rbar = norm(sumx)/n;
% muMLE = sumx/norm(sumx);
kappaMLE = (rbar*d - rbar^3) / (1-rbar^2);
[centers, rs] = getSubSphere(data,0);
radialdistances = acos(centers'*data);
xi_sample = mean(radialdistances)/std(radialdistances);

% 2. Now generate boostrap samples from \mu = muMLE , \kappa = kappaMLE
% (But WLOG, \mu can be assumed to be e_1)
xi_vec = zeros(R,1);
%tstart = tic;

if (g_use_parallel)
    parfor r = 1:R
        rdata = randvonMisesFisherm(d,n,kappaMLE);
        [centers, rs] = getSubSphere(rdata,0);
        radialdistances = acos(centers'*rdata);
        xi_vec(r) = mean(radialdistances)/std(radialdistances);
    end
else
    for r = 1:R
        rdata = randvonMisesFisherm(d,n,kappaMLE);
        [centers, rs] = getSubSphere(rdata,0);
        radialdistances = acos(centers'*rdata);
        xi_vec(r) = mean(radialdistances)/std(radialdistances);
    end
end

%telapsed = toc(tstart);
%fprintf('\n');
%fprintf('--------------------------\n');
%fprintf('vMftest:elapsed time = %f\n',telapsed);
%if g_use_parallel; disp(' parfor used'); end;
%fprintf('--------------------------\n\n');

% Now see how many boostrap samples are > xi_sample
pvalue = mean(xi_vec > xi_sample);

% % Bonus: Let's see the null distribution
% clf;
% [kde,xgrid,mker] = kdeSM(xi_vec,struct('ibdryadj',0,'iplot',1));
% xlabel(['kappa ' num2str(kappa) ', ratios estimated']);title(['S^d, d = ' num2str(d) ', n = ' num2str(n)]);
% ylabel(['pvalue' num2str(pvalue)]);
% hold on; 
% plot([xi_sample xi_sample],[0 max(kde)],'-k','Linewidth',0.5);

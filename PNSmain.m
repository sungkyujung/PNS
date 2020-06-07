function [resmat, PNS]= PNSmain(data,itype,alpha,R,parflag)
% PNSMAIN Analysis of Principal Nested Spheres for data on hyperspheres
% [resmat PNS]= PNSmain(data)
% [resmat PNS]= PNSmain(data,itype)
% [resmat PNS]= PNSmain(data,itype,alpha)
% [resmat PNS]= PNSmain(data,itype,alpha,R)
% [resmat PNS]= PNSmain(data,itype,alpha,R,parflag)
% 
% 
% Input:
%   data     (d+1) x n data matrix where each column is a unit vector.
%
%   itype    0  'seq.test' : (default) ordinary Principal Nested Sphere
%                               with sequential tests.
%            1  'small'    : Principal Nested SMALL Sphere
%            2  'great'    : Principal Nested GREAT Sphere (radius pi/2)
%            3  'BIC'      : BIC rule
%            4  'RS'       : Rotational symmetry test 
%            6  'S1'       : Use the small sphere distribution of the first
%                               kind (S1) and related likelihood ratio test
%                               in sequential tests.
%                            Requires fitS1.m and fitS1givennu.m functions.
%            7  'FN'  :  Folded normal test by Eltzner et al. only. 
%            8  'FN + LRT'  : Goodness-of-fit test + Folded normal test
%            9  'Kurtosis': Goodness-of-fit test + Kurtosis test of modality
%            10 'BIC+Kurt': BIC rule + Kurtosis test of modality
%            11 'BIC+FN  ': BIC rule + Folded normal test
%
%  alpha     0.05 (default) : size of Type I error allowed for each test,
%            could be any number between 0 and 1.
%  R         100 (default) : number of bootsrap samples to be evaluated for
%            the sequential test.
%
%  parflag  false (default) : If parflag is set to "true", parfor is
%                             enabled in vMFTest.m.
%
% Output:
%   resmat  The commensurate residual matrix (X_PNS). Each entry in row k
%           works like the kth principal c[omponent score.
%
%   PNS     Matlab structure array with following fields:
%       mean     : location of the PNSmean
%       radii    : size (radius) of PNS
%       orthoaxis: orthogonal axis 'v_i' of subspheres
%       dist     : distance 'r_i' of subspheres
%       pvalues  : p-values of the corresponding tests
%       gsphere  : indicator for great (=1) or small (=0) sphere fitting
%       itype    : type of methods for fitting subspheres
%
% Last updated October 2019
% Sungkyu Jung,Takayuki Iwamoto and Byungwon Kim%
%
% See also PNSshape

%% Handling input parameters 
global g_use_parallel;

if nargin == 1;
    itype = 0; %  'seq.test'; % default: sequential test
    alpha = 0.05;
    R = 100;
end

if nargin == 2;
    alpha = 0.05;
    R = 100;
end

if nargin == 3;
    R = 100;
end


if nargin < 5;
    if g_use_parallel % if g_use_parallel == true
    else              % all other cases;
        g_use_parallel = false;
    end
end

if nargin == 5;
    % if parflage is read, override g_use_parallel = parflag
    if parflag % if parflag == true
        % try open matlabpool if not open yet;
        try
            %             if (parpool('size') == 0)
            %                 parpool;
            %             end
            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(poolobj)
                parpool;
            end
        catch
            parflag  = false;
        end
    end
    g_use_parallel = parflag; % true or false
end



%% Data matrix handling 

[k, n] = size(data); % data on (k-1)-sphere in Real^k
[uu, ll]=svd(data,'econ');
maxd = find(diag(ll) < 1e-15,1);
if isempty(maxd) || k > n
    maxd = min(k,n)+1;
end
nullspdim = k - maxd + 1; % dimension of subspace that contains no data

d = k - 1; % intrinsic dimension of the sphere
disp(['Message from PNSmain.m; dataset is on ' num2str(d) '-sphere.']);

dm = maxd - 2; % intrinsic dimension of the smallest nested sphere that
% contains variation. I.e. this is the dimension to which
% can be trivially reduced.

if k == 3
    nullspdim = 0;
    dm = 2;
end


resmat = zeros(dm,n); % d dimensional residual matrix
% there will be dm-1 subspheres
orthaxis = cell(dm-1,1);
dist = zeros(dm-1,1);
%pvalues = zeros(dm-1,2);
pvalues = NaN(dm-1,2);
gsphere = NaN(dm-1,1);
iso = NaN(dm-1,1);

if nullspdim > 0
    disp([' .. found null space of dimension ' num2str(nullspdim)...
        ',to be trivially reduced.']);
    disp([' .. then narrow down to ' num2str(dm) '-sphere']);
    % (HDLSS case) fit nested great spheres for dimension reduction
    % where no residual is present.
    currentSphere = uu(:,1:(dm+1))'*data ;
else
    disp([' .. Check that the following holds: ' num2str(d) ' = ' ...
        num2str(dm)]);
    currentSphere = data;
end



if itype == 0;  % 'seq.test'
    disp([' .. sequential tests with significance level ' num2str(alpha)])
    isIsotropic = false; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0); % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs;
                [centerg, rg] = getSubSphere(currentSphere,1); % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg;
                
                % Chi-squared statistic for a likelihood test
                pval1 = LRTpval(resGREAT,resSMALL,n);
                
                pvalues(i,1) = pval1;
                if pval1 > alpha
                    center = centerg; r = rg;
                    pvalues(i,2) = NaN;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ', p(LRT) = ' num2str(pval1)]);
                else
                    pval2 = vMFtest(currentSphere,R);
                    pvalues(i,2) = pval2;
                    if pval2 > alpha
                        center = centerg; r = rg;
                        gsphere(i) = 1;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'GREAT sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(vMF) = ' num2str(pval2)]);
                        isIsotropic = true;
                    else
                        center = centers; r = rs;
                        gsphere(i) = 0;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'SMALL sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(vMF) = ' num2str(pval2)]);
                    end
                end
                
            case true
                [center, r] = getSubSphere(currentSphere,1);
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by testing vMF distn']);
                pvalues(i,1) = NaN;
                pvalues(i,2) = NaN;
        end
        iso(i) = isIsotropic;
        res = acos(center'*currentSphere)-r;
        % save subsphere parameters
        orthaxis{i} = center; dist(i) = r;
        % save residuals
        resmat(i,:) = res;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]);
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1);
    end
    
elseif itype == 3;  % 'BIC'
    disp(' .. with BIC')
    for i=1:(dm-1) %i=nullspdim+1:(d-1)
        [centers, rs] = getSubSphere(currentSphere,0); % small sphere fit
        resSMALL = acos(centers'*currentSphere)-rs;
        [centerg, rg] = getSubSphere(currentSphere,1); % great sphere fit
        resGREAT = acos(centerg'*currentSphere)-rg;
        
        % BIC
        BICsmall = n*log(mean(resSMALL.^2)) + (dm-i+1+1)*log(n);
        BICgreat = n*log(mean(resGREAT.^2)) + (dm-i+1)*log(n);
        disp(['BICsm: ' num2str(BICsmall) ', BICgr: ' num2str(BICgreat)])
        if BICsmall > BICgreat
            center = centerg; r = rg;
            gsphere(i) = 1;
            disp([num2str(dm-i+1) '-sphere to '...
                num2str(dm-i) '-sphere, by '...
                'GREAT sphere' ', BIC ' ]);
        else    center = centers; r = rs;
	    gsphere(i) = 0;
            disp([num2str(dm-i+1) '-sphere to '...
                num2str(dm-i) '-sphere, by '...
                'SMALL sphere' ...
                ', BIC' ]);
        end
        res = acos(center'*currentSphere)-r;
        % save subsphere parameters
        orthaxis{i} = center; dist(i) = r;
        % save residuals
        resmat(i,:) = res;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1);
    end
    
    
elseif itype == 4;  % 'Tangent space rotational symmetry'
    disp([' .. sequential tests (Rotational Symmetry) with significance level ' num2str(alpha)])
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0) ; % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs ;
                [centerg, rg] = getSubSphere(currentSphere,1) ; % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg ;

                % Chi-squared statistic for a likelihood test
                pval1 = LRTpval(resGREAT,resSMALL,n) ;

                pvalues(i,1) = pval1 ;
                if pval1 > alpha
                    center = centerg ; r = rg ;
                    pvalues(i,2) = NaN ;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ', p(LRT) = ' num2str(pval1)]) ;
                else
                    pval2 = TangentNormaltest(currentSphere) ;
                    pvalues(i,2) = pval2 ;
                    if pval2 > alpha
                        center = centerg ; r = rg ;
                        gsphere(i) = 1;
                        disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by '...
                                'GREAT sphere' ', p(LRT) = ' num2str(pval1)...
                                ', p(TN) = ' num2str(pval2)]);
                            isIsotropic = true ;
                    else
                        center = centers ; r = rs ;
                        gsphere(i) = 0;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'SMALL sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(TN) = ' num2str(pval2)]) ;
                    end
                end

            case true
                [center, r] = getSubSphere(currentSphere,1) ;
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by testing tangent normal']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        iso(i) = isIsotropic;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end


elseif itype == 6;  % 'S1'
    disp([' .. sequential tests by S1 with significance level ' num2str(alpha)])
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                % testing for great circle (nu = 0)
                S1_Great = fitS1givennu(currentSphere, 0) ;
                S1_Small = fitS1(currentSphere) ;

                tstat = max(2 * (S1_Great.nll - S1_Small.nll), 0) ;
                pval1 = chi2cdf(tstat, 1, 'upper') ;
                
                pvalues(i,1) = pval1 ;
                if pval1 > alpha
                    [center, r] = getSubSphere(currentSphere,1) ; % great sphere fit
                    pvalues(i,2) = NaN ;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ', p(S1LRT) = ' num2str(pval1)]) ;
                else
                    % testing for isotropy (kappa_0 = 0)
                    Xbar = mean(currentSphere, 2) ;
                    Rbar = norm(Xbar) ;
                    mhat = Xbar / Rbar ;
                    khat = Rbar * (3 - Rbar^2) / (1 - Rbar^2) ; % Banerjee et al (2005)
                    nllvMF = -n*log(khat) + n*log(4*pi*sinh(khat)) - n*khat*mhat'*Xbar ; % negative log likelihood

                    tstat = max(2 * (nllvMF - S1_Small.nll), 0) ;
                    pval2 = chi2cdf(tstat, (dm-i+2), 'upper') ;
    
                    pvalues(i,2) = pval2 ;
                    if pval2 > alpha
                        [center, r] = getSubSphere(currentSphere,1) ; % great sphere fit
                        gsphere(i) = 1;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'GREAT sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(S1vMF) = ' num2str(pval2)]);
                        isIsotropic = true ;
                    else
                        [center, r] = getSubSphere(currentSphere,0) ; % small sphere fit
                        gsphere(i) = 0;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'SMALL sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(S1vMF) = ' num2str(pval2)]) ;
                    end
                end

            case true
                [center, r] = getSubSphere(currentSphere,1) ;
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by testing S1']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        iso(i) = isIsotropic;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end

elseif itype == 7;  % 'Eltzner'
    disp([' .. sequential tests by Eltzner et al. (2015+) with significance level ' num2str(alpha)])
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [center, r] = getSubSphere(currentSphere,0) ;
                radii = acos(center'*currentSphere) ;
                Mean = mean(radii) ;
                Std = std(radii) ;
                f = @(r, rho, sigma) (exp(-0.5*(r / sigma - rho).^2) ...
                                     + exp(-0.5*(r / sigma + rho).^2)) ...
                                     .* (sin(r).^(dm-i+2-1)) ;
                normalization = @(rho, sigma) integral(@(r)f(r,rho,sigma), 0, pi) ;
                likelihood = @(param) n * log(normalization(param(1), param(2))) ...
                                           - (dm-i+2-1) * sum(log(sin(radii))) ...
                                           + 0.5 * sum(radii.^2 / param(2)^2 + param(1)^2) ...
                                           - sum(log(cosh(radii*param(1)/param(2)))) ;
                [mle,fval,exitflag,output] = fminsearchbnd(likelihood, [Mean/Std;Std], [0;0], [inf;inf]) ;
                rhohat = mle(1) ;
                sigmahat = mle(2) ;
                if rhohat > 1 ;
                    likelihood_null = @(sigma) likelihood([1; sigma]) ;
                    [mle_null,fval,exitflag,output] = fminsearchbnd(likelihood_null, Std, 0, inf) ;
                    sigmahat_null = mle_null(1) ;

                    pval = chi2cdf(-2*(likelihood([rhohat; sigmahat]) - likelihood_null(sigmahat_null)), 1, 'upper') ;
                    pvalues(i,1) = pval ;
                    pvalues(i,2) = NaN ;
                    if pval > alpha ;
                        [center, r] = getSubSphere(currentSphere,1) ; % great sphere fit
                        gsphere(i) = 1;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'GREAT sphere' ', p(LRT(rho=1)) = ' num2str(pval)]) ;
                    else
                        [center, r] = getSubSphere(currentSphere,0) ; % small sphere fit
                        gsphere(i) = 0;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'SMALL sphere' ', p(LRT(rho=1)) = ' num2str(pval)]) ;
                    end
                else
                    pvalues(i,1) = NaN ;
                    pvalues(i,2) = NaN ;
                    [center, r] = getSubSphere(currentSphere,1) ; % great sphere fit
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ', estimated rho is less than 1 (Isotropic)']) ;
                    isIsotropic = true ;
                end

            case true
                [center, r] = getSubSphere(currentSphere,1) ;
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by test of Eltzner et al (2015+)']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end
    
elseif itype == 8;  % 'Eltzner + LRT'
    disp([' .. sequential tests by Eltzner et al. (2015+)' ...
          ' + great sphere LRT with significance level ' num2str(alpha)])
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0) ; % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs ;
                [centerg, rg] = getSubSphere(currentSphere,1) ; % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg ;

                % Chi-squared statistic for a likelihood test
                pval1 = LRTpval(resGREAT,resSMALL,n) ;

                pvalues(i,1) = pval1 ;
                if pval1 > alpha
                    center = centerg ; r = rg; % great sphere fit
                    pvalues(i,2) = NaN;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ', p(LRT(rho*sigma=pi/2)) = ' num2str(pval1)]);
                else
		    center = centers ; r = rs;
		    radii = acos(center'*currentSphere) ;
		    Mean = mean(radii) ;
		    Std = std(radii) ;
		    f = @(r, rho, sigma) (exp(-0.5*(r / sigma - rho).^2) ...
					+ exp(-0.5*(r / sigma + rho).^2)) ...
					.* (sin(r).^(dm-i+2-1)) ;
		    normalization = @(rho, sigma) integral(@(r)f(r,rho,sigma), 0, pi) ;
		    likelihood = @(param) n * log(normalization(param(1), param(2))) ...
					    - (dm-i+2-1) * sum(log(sin(radii))) ...
					    + 0.5 * sum(radii.^2 / param(2)^2 + param(1)^2) ...
					    - sum(log(cosh(radii*param(1)/param(2)))) ;
		    [mle,fval,exitflag,output] = fminsearchbnd(likelihood, [Mean/Std;Std], [0;0], [inf;inf]) ;
		    rhohat = mle(1) ;
		    sigmahat = mle(2) ;

		    
                    if rhohat > 1 ;
                        likelihood_null2 = @(sigma) likelihood([1; sigma]) ;
                        [mle_null,fval,exitflag,output] = fminsearchbnd(likelihood_null2, Std, 0, inf) ;
                        sigmahat_null = mle_null(1) ;

                        pval2 = chi2cdf(-2*(likelihood([rhohat; sigmahat]) - likelihood_null2(sigmahat_null)), 1, 'upper') ;
                        pvalues(i,2) = pval2 ;
                        if pval2 > alpha ;
                            center = centerg ; r = rg; % great sphere fit
                            gsphere(i) = 1;
                            disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by ' 'GREAT sphere' ...
                                ', p(LRT(rho*sigma=pi/2)) = ' num2str(pval1) ...
                                ', p(LRT(rho=1)) = ' num2str(pval2)]) ;
                        else
                            center = centers ; r = rs; % small sphere fit
                            gsphere(i) = 0;
                            disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by ' 'SMALL sphere' ...
                                ', p(LRT(rho*sigma=pi/2)) = ' num2str(pval1) ...
                                ', p(LRT(rho=1)) = ' num2str(pval2)]) ;
                        end
                    else
                        pvalues(i,2) = NaN ;
                        center = centerg ; r = rg; % great sphere fit
                        gsphere(i) = 1;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by ' 'GREAT sphere' ...
                            ', p(LRT(rho*sigma=pi/2)) = ' num2str(pval1) ...
                            ', estimated rho is less than 1 (Isotropic)']) ;
                        isIsotropic = true ;
                    end
                end

            case true
                [center, r] = getSubSphere(currentSphere,1) ;
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by test of Eltzner et al (2015+)']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        iso(i) = isIsotropic;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end

elseif itype == 9;  % 'GF + Modality test'
    disp([' .. sequential tests (Modality test) with significance level ' num2str(alpha)])
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0) ; % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs ;
                [centerg, rg] = getSubSphere(currentSphere,1) ; % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg ;

                % Chi-squared statistic for a likelihood test
                pval1 = LRTpval(resGREAT,resSMALL,n) ;

                pvalues(i,1) = pval1 ;
                if pval1 > alpha
                    center = centerg ; r = rg ;
                    pvalues(i,2) = NaN ;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ', p(LRT) = ' num2str(pval1)]) ;
                else
                    
                    % %%% kurtosis test routine %%% %
                    X = LogNPd(rotMat(geodmeanSk(centers)) * currentSphere) ;
                    % Note that the tangential point is the center of the small circle
                    [d, n] = size(X) ;
                    normX2 = sum(X.^2);
                    kurtosis = sum( normX2.^2 ) / n / ( sum( normX2 ) / (d * (n-1)) )^2;
                    M_kurt = d * (d+2)^2 / (d+4) ;
                    V_kurt = (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8)) ;
                    pval2 = normcdf((kurtosis - M_kurt) / sqrt(V_kurt)) ; 
                    
                    % %%% kurtosis test routine - end %%% % 
                    pvalues(i,2) = pval2 ;
                    if pval2 > alpha
                        center = centerg ; r = rg ;
                        gsphere(i) = 1;   
                        disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by '...
                                'GREAT sphere' ', p(LRT) = ' num2str(pval1)...
                                ', p(TN) = ' num2str(pval2)]);
                            isIsotropic = true ;
                    else
                        center = centers ; r = rs ;
                        gsphere(i) = 0;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'SMALL sphere' ', p(LRT) = ' num2str(pval1)...
                            ', p(TN) = ' num2str(pval2)]) ;
                    end
                end

            case true
                [center, r] = getSubSphere(currentSphere,1);
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by testing tangent normal']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        iso(i) = isIsotropic;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end
 
elseif itype == 10;  % 'BIC + Modality test'
    disp([' .. sequential decision by BIC + Modality test']) 
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0) ; % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs ;
                [centerg, rg] = getSubSphere(currentSphere,1) ; % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg ;
                
                % BIC
                BICsmall = n*log(mean(resSMALL.^2)) + (d-i+1+1)*log(n);
                BICgreat = n*log(mean(resGREAT.^2)) + (d-i+1)*log(n);
                
                pvalues(i,1) = BICsmall - BICgreat ;
                if BICsmall > BICgreat
                    center = centerg ; r = rg ;
                    pvalues(i,2) = NaN ;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ...
                        'since BICsm: ' num2str(BICsmall) ' >  BICgr: ' num2str(BICgreat)]) ;
                else
                    
                    % %%% kurtosis test routine %%% %
                    X = LogNPd(rotMat(geodmeanSk(centers)) * currentSphere) ;
                    % Note that the tangential point is the center of the small circle
                    [d, n] = size(X) ;
                    normX2 = sum(X.^2);
                    kurtosis = sum( normX2.^2 ) / n / ( sum( normX2 ) / (d * (n-1)) )^2;
                    M_kurt = d * (d+2)^2 / (d+4) ;
                    V_kurt = (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8)) ;
                    pval2 = normcdf((kurtosis - M_kurt) / sqrt(V_kurt)) ; 
                    
                    % %%% kurtosis test routine - end %%% % 
                    pvalues(i,2) = pval2 ;
                    if pval2 > alpha
                        center = centerg ; r = rg ;
                        gsphere(i) = 1;   
                        disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by '...
                                'GREAT sphere' ',BICsm: ' num2str(BICsmall) ' < BICgr: ' num2str(BICgreat) ...
                                'but p(Modality) = ' num2str(pval2)]);
                            isIsotropic = true ;
                    else
                        center = centers ; r = rs ;
                        gsphere(i) = 0;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by '...
                            'SMALL sphere' ',BICsm: ' num2str(BICsmall) '< BICgr: ' num2str(BICgreat) ...
                            ', p(Modality) = ' num2str(pval2)]) ;
                    end
                end

            case true
                [center, r] = getSubSphere(currentSphere,1);
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by Kurtosis test of modality']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        iso(i) = isIsotropic;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end

elseif itype == 11;  % 'Eltzner + LRT'
    disp([' .. sequential decision by BIC + Folded normal test']) 
    isIsotropic = false ; % to be used as a test procedure
    for i= 1:(dm-1)
        switch isIsotropic
            case false
                [centers, rs] = getSubSphere(currentSphere,0) ; % small sphere fit
                resSMALL = acos(centers'*currentSphere)-rs ;
                [centerg, rg] = getSubSphere(currentSphere,1) ; % great sphere fit
                resGREAT = acos(centerg'*currentSphere)-rg ;

                % BIC
                BICsmall = n*log(mean(resSMALL.^2)) + (d-i+1+1)*log(n);
                BICgreat = n*log(mean(resGREAT.^2)) + (d-i+1)*log(n);
                

                pvalues(i,1) = BICsmall - BICgreat ;
                if BICsmall > BICgreat
                    center = centerg ; r = rg ;
                    pvalues(i,2) = NaN ;
                    gsphere(i) = 1;
                    disp([num2str(dm-i+1) '-sphere to '...
                        num2str(dm-i) '-sphere, by '...
                        'GREAT sphere' ...
                        'since BICsm: ' num2str(BICsmall) ' >  BICgr: ' num2str(BICgreat)]) ;
                else
		    center = centers ; r = rs;
		    radii = acos(center'*currentSphere) ;
		    Mean = mean(radii) ;
		    Std = std(radii) ;
		    f = @(r, rho, sigma) (exp(-0.5*(r / sigma - rho).^2) ...
					+ exp(-0.5*(r / sigma + rho).^2)) ...
					.* (sin(r).^(dm-i+2-1)) ;
		    normalization = @(rho, sigma) integral(@(r)f(r,rho,sigma), 0, pi) ;
		    likelihood = @(param) n * log(normalization(param(1), param(2))) ...
					    - (dm-i+2-1) * sum(log(sin(radii))) ...
					    + 0.5 * sum(radii.^2 / param(2)^2 + param(1)^2) ...
					    - sum(log(cosh(radii*param(1)/param(2)))) ;
		    [mle,fval,exitflag,output] = fminsearchbnd(likelihood, [Mean/Std;Std], [0;0], [inf;inf]) ;
		    rhohat = mle(1) ;
		    sigmahat = mle(2) ;

		    
                    if rhohat > 1 ;
                        likelihood_null2 = @(sigma) likelihood([1; sigma]) ;
                        [mle_null,fval,exitflag,output] = fminsearchbnd(likelihood_null2, Std, 0, inf) ;
                        sigmahat_null = mle_null(1) ;

                        pval2 = chi2cdf(-2*(likelihood([rhohat; sigmahat]) - likelihood_null2(sigmahat_null)), 1, 'upper') ;
                        pvalues(i,2) = pval2 ;
                        if pval2 > alpha ;
                            center = centerg ; r = rg; % great sphere fit
                            gsphere(i) = 1;
                            disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by ' 'GREAT sphere' ...
                                ',BICsm: ' num2str(BICsmall) ' < BICgr: ' num2str(BICgreat) ...
                                ', p(LRT(rho=1)) = ' num2str(pval2)]) ;
                        else
                            center = centers ; r = rs; % small sphere fit
                            gsphere(i) = 0;
                            disp([num2str(dm-i+1) '-sphere to '...
                                num2str(dm-i) '-sphere, by ' 'SMALL sphere' ...
                                ',BICsm: ' num2str(BICsmall) ' < BICgr: ' num2str(BICgreat) ...
                                ', p(LRT(rho=1)) = ' num2str(pval2)]) ;
                        end
                    else
                        pvalues(i,2) = NaN ;
                        center = centerg ; r = rg; % great sphere fit
                        gsphere(i) = 1;
                        disp([num2str(dm-i+1) '-sphere to '...
                            num2str(dm-i) '-sphere, by ' 'GREAT sphere' ...
                            ',BICsm: ' num2str(BICsmall) ' < BICgr: ' num2str(BICgreat) ...
                            ', estimated rho is less than 1 (Isotropic)']) ;
                        isIsotropic = true ;
                    end
                end

            case true
                [center, r] = getSubSphere(currentSphere,1) ;
                gsphere(i) = 1;
                disp([num2str(dm-i+1) '-sphere to '...
                    num2str(dm-i) '-sphere, by '...
                    'GREAT sphere, restricted by Folded normal test of Eltzner et al']) ;
                pvalues(i,1) = NaN ;
                pvalues(i,2) = NaN ;
        end
        res = acos(center'*currentSphere)-r ;
        iso(i) = isIsotropic;
        % save subsphere parameters
        orthaxis{i} = center ; dist(i) = r ;
        % save residuals
        resmat(i,:) = res ;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere ;
        %disp(['Check that ' num2str(size(NestedSphere)) ' = ' num2str(dm+2-i)]) ;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1) ;
    end

    
elseif itype == 1|| itype ==2 % strictly 'small' or 'great' spheres
    pvalues = NaN;
    for i= 1:(dm-1)
        % estimate the best fitting subsphere
        %  with small sphere if itype = 1
        %  with  great sphere if itype = 2
        [center, r] = getSubSphere(currentSphere,itype-1);
        res = acos(center'*currentSphere)-r;
        % save subsphere parameters
        orthaxis{i} = center;
        dist(i) = r;
        % save residuals
        resmat(i,:) = res;
        % projection to subsphere and transformation to isomorphic
        % sphere
        NestedSphere = rotMat(center)*currentSphere;
        currentSphere = NestedSphere(1:(dm+1-i),:)./...
            repmat(sqrt(1-NestedSphere(dm+2-i,:).^2),dm+1-i,1);
        gsphere(i) = itype-1;
    end
else
    disp('!!! Error from PNSmain.m:');
    disp('!!! itype must be 0 (seq.test), 1 (small), 2 (great) or up to 10.');
    disp('!!!   Terminating execution     ') ;
    return ;
end

% currentSphere has intrinsic dimension 1
% compute PNSmean and deviations.

% parametrize 1-sphere to angles
if nullspdim+1-(dm-1) <= 0
    S1toRadian = atan2(currentSphere(2,:),currentSphere(1,:));
    % Geodesic mean of angles
    meantheta = geodmeanS1(S1toRadian');
    orthaxis{dm} = meantheta;
    % save deviations from PNSmean
    resmat(dm,:) =mod(S1toRadian - meantheta + pi,2*pi)-pi;
end

%% when there is no variation. 3/18/2013.
if nullspdim+1-(dm-1) > 0
    S1toRadian = atan2(currentSphere(2,:),currentSphere(1,:));
    % Geodesic mean of angles
    meantheta = geodmeanS1(S1toRadian');
    orthaxis{dm} = meantheta;
    % save deviations from PNSmean
    resmat(dm,:) =mod(S1toRadian - meantheta + pi,2*pi)-pi;
end
%%

radii =1;
for i = 1:(dm-1)
    radii = [radii; prod( sin( dist(1:i) ) )];
end
resmat = flipud(repmat(radii,1,n).*resmat); % scaled residuals

PNS.radii = radii;          % size (radius) of nested spheres from largest to smallest
PNS.orthaxis = orthaxis;    % orthogonal axis of (d-1) subspheres and the anglemean for PNSmean
PNS.dist = dist;            % distances for (d-1) subspheres
PNS.pvalues = pvalues;      % d-1 pvalues from sequential tests
PNS.gsphere = gsphere;      % d-1 indicator for great (1) or small (0) sphere

if nullspdim > 0;
    PNS.basisu = uu(:,1:(dm+1));
else
    PNS.basisu = [];
end

PNS.mean = PNSe2s(zeros(dm,1),PNS); % PNSmean of the data

switch itype
    case 0
        PNS.itype = 'seq.test';
    case 1
        PNS.itype = 'small';
    case 2
        PNS.itype = 'great';
end

end

function pval = LRTpval(resGREAT,resSMALL,n)
chi2 = max(n*log(sum(resGREAT.^2)/sum(resSMALL.^2)),0);
pval = 1-chi2cdf(chi2,1); % p-value of the likelihood test
end

function pval = TangentNormaltest(currentSphere)
% Testing H0: Sigma = c*I where c>0 on the tangent space of geodesic mean

dataontangent = LogNPd(rotMat(geodmeanSk(currentSphere)) * currentSphere) ;

[D N] = size(dataontangent) ;
Mest = mean(dataontangent, 2) ; % Under both H0 and H1
SestH1 = cov(dataontangent') + 1e-5*eye(D) ; % Under H1
 SestH0 = sum(diag((dataontangent - repmat(Mest,1,N))' * ...
                   (dataontangent - repmat(Mest,1,N)))) ...
          / (D*(N-1)) * eye(D) ;
 Teststat = -2 * ( sum(log(mvnpdf(dataontangent', Mest', SestH0))) ...
                   - sum(log(mvnpdf(dataontangent', Mest', SestH1))) ) ;
pval = 1 - chi2cdf(Teststat, D*(D+1)/2-1) ;
end 

function [kurtosis, pval] = TestBallUniform(X)
% TestBallUniform(X) returns the p-value from the test of uni-modality
% using the asymptotic normal distribution of sample multivariate kurtosis
% obtained from assuming the underlying distribution to be a ball uniform
% distribution.
%
% Input - X: (d x n) data matrix, d is dimension, n is sample size
% Output - kurtosis: the measure of multivariate kurtosis
%                    introduced by Mardia (1970)
%          pval: the p-value calculated for the test of uni-modality
%
% Requires MardiaKurtosis.m

[d, n] = size(X) ;

%% Calculate the sample kurtosis
[kurtosis] = MardiaKurtosis(X, 1) ;
% we use the option '1' for this function to estimate the covariance matrix
% with a restriction, \sigma^2*I_d. <- We assume the distribution is
% isotropic (rotationally symmetric with respect to mean)
M_kurt = d * (d+2)^2 / (d+4) ;
V_kurt = (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8)) ;
pval = normcdf((kurtosis - M_kurt) / sqrt(V_kurt)) ;
end

function [Kurtosis] = MardiaKurtosis(X, option)
% This function calculates the measure of multivariate kurtosis introduced
% by Mardia (1970)
%
% [kurtosis] = MardiaKurtosis(X, option)
% X: data matrix (d x n) 
% option: 0 (default) - uses the common covariance matrix
%         1 - Assume the distribution is isotropic
%             With this option, the covariance matrix is estimated by
%             assuming the form of covariance matrix as (scalar x identity)
%
% Byungwon Kim, Sep 2017
if nargin == 1
    option = 0 ;
end

[D, N] = size(X) ;
M = mean(X, 2) ;
if option == 0 ;
    S = cov(X') + 1e-5*eye(D) ;
else
    S = sum(diag((X - repmat(M,1,N))' * (X - repmat(M,1,N)))) / (D*(N-1)) ...
        * eye(D) ;
end
Kurtosis = sum(diag((X - repmat(M,1,N))' / S * (X - repmat(M,1,N))).^2) ...
           / N ;
end
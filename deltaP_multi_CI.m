% -----------------------------------------------------------------------------------------------------------
% function [deltaPs, rawPercentiles] = deltaP_multi_CI(inputData, nResamples = 1000, alpha = 0.05, method = [3, 1], useMedian = true)
% -----------------------------------------------------------------------------------------------------------
% Confidence intervals for deltaP (version for multiple traits);
% inputData must be in the stacked format (matrix with trait values in
% the first column, population labels in the second column, and trait
% labels in the third column); if the third column is missing, 
% a single trait is assumed;
% nResamples is the number of bootstrap resamples;
% 1 - alpha is the confidence level;
% useMedian indicates whether deltaP is calculated based on population
% medians (default) or means;
% method is a vector indicates the methods used for constructing CIs; 
% the first entry gives the method for the univariate deltaP (single-trait comparisons),
% and the second entry gives the method for the multivariate deltaP (Euclidean 
% distance in percentile space; 
% method == 1: percentile method
% method == 2: BC method (bias-corrected)
% method == 3: BCa method (bias-corrected and accellerated
% defaults are BCa (3) for univariate deltaP and percentile (1) for multivariate deltaP;
% The output is a matrix with rows corresponding to the individual
% comparisons:
% column 1: label of population 1
% column 2: label of population 2
% column 3: trait label (-1 for multivariate comparison at the end of the ouptut)
% column 4: point estimate for deltaP
% column 5: the bias-corrected point estimate: 2 DeltaP minus mean(bootstrap-distribution)
%	    this correction is only applied for the multivariate comparison
% column 6: lower confidence limit (for the chosen method)
% column 7: upper confidence limit
% column 8: p-value from permutation test
function [y, P] = deltaP_multi_CI(inputData, nResamples, alpha, method, useMedian)
    if nargin < 5
        useMedian = true;
    end
    if nargin < 4
        method = [3, 1];
    end
    if nargin < 3
        alpha = 0.05;
    end
    if nargin < 2
        nResamples = 1000;
    end
    if size(inputData, 2) == 2
        inputData(:,3) = 1;
    end
    if size(inputData, 2) ~= 3
        error('Input data must be a n x 2 or n x 3 matrix.\n');
    end
    disp(' ');
    disp('...Working...');
    disp('...Using bootstrapping method to calculate confidence intervals for delta-P values...');
    disp(' ');
    [temp,P] = deltaP_multi(inputData, useMedian);
    sampleValue = temp(:,4)'; 
    temp = deltaP_multi_permutation(inputData, nResamples, useMedian);
    pValue = temp(:,5)';
    data = inputData(:,1);
    popVector = inputData(:,2);
    popLabels = unique(popVector);
    nPops = length(popLabels);
    traitVector = inputData(:,3);
    traitLabels = unique(traitVector);
    nTraits = length(traitLabels);
    % create bootstrap distribution
    samplesizes = zeros(nTraits,nPops);
    for i = 1:nTraits
	for j = 1:nPops
	    sample{i,j} = data( traitVector==traitLabels(i) & popVector==popLabels(j) );
	    samplesizes(i,j) = length(sample{i,j});
	end
    end
    bootdist = zeros(nResamples, (nTraits+1)*nPops*(nPops-1)/2);
    for k = 1:nResamples
	resample = [];
	for i = 1:nTraits
	    for j = 1:nPops
		boot = ceil(samplesizes(i,j) * rand(1, samplesizes(i,j)));
		clear s;
		s = sample{i,j}(boot);
		s(:,2) = popLabels(j);
		s(:,3) = traitLabels(i);
		resample = [resample; s];
	    end
	end
	[temp,~] = deltaP_multi(resample, useMedian);
	bootdist(k, :) = temp(:,4)';
    end
    % create jackknife distributions (BCa method only)
    % single-trait comparisons
    if method(1) == 3
	for j = 1:nTraits
	    tSample = inputData(traitVector == traitLabels(j), 1:2);
	    jDist = zeros(length(tSample), nPops*(nPops-1)/2);
	    for i = 1:length(tSample)
		jSample = tSample;
		jSample(i,:) = [];
		temp = deltaP_single(jSample, useMedian);
		jDist(i,:) = temp(:,3)';
	    end
	    jackdist{j} = jDist;
	end
    end
    % multivariate comparisons (Euclidean distances)
    if method(2) == 3
	nComp = nPops*(nPops-1)/2;
	jDist = zeros(length(data), nComp);
	for i = 1:length(data)
	    jSample = inputData;
	    jSample(i,:) = [];
	    [temp,~] = deltaP_multi(jSample, useMedian);
	    jDist(i,:) = temp(end-nComp+1:end,4)';
	end
	jackdist{nTraits+1} = jDist;
    end
    % calculate BCa CIs 
    y = zeros((nTraits+1)*nPops*(nPops-1)/2, 8);
    h = 0;
    for j = 1:nTraits+1
	h1 = 0;
	for i1 = 1:nPops
	    for i2 = (i1+1):nPops
		h = h + 1;
		h1 = h1 + 1;
		% acceleration (BCa method only)
		if ( (j <= nTraits) && (method(1) == 3) ) || ( (j > nTraits) && (method(2) == 3) )
		    meanJack = mean(jackdist{j}(:,h1));
		    if (var(jackdist{j}(:,h1)) == 0)
			a = 0;
		    else
			a = (sum((meanJack-jackdist{j}(:,h1)).^3))/(6*sum((meanJack-jackdist{j}(:,h1)).^2).^(3/2)); 
		    end
		else
		    a = 0;
		end
		% bias correction (BC and BCa methods only)
		if ( (j <= nTraits) && (method(1) >= 2) ) || ( (j > nTraits) && (method(2) >= 2) )
		    z0 = norminv( ( sum(bootdist(:,h) < sampleValue(h)) + sum(bootdist(:,h) == sampleValue(h))/2. ) / nResamples); 
		else
		    z0 = 0;
		end
		% calculate confidence limits
		if z0 == Inf 
		    CILowBCa = 1;
		    CIUpBCa = 1;
		elseif z0 == -Inf
		    CILowBCa = 0;
		    CIUpBCa = 0;
		else
		    CILowBCa = normcdf(z0 + (z0+norminv(alpha/2))/(1-a*(z0+norminv(alpha/2))));
		    CIUpBCa  = normcdf(z0 + (z0+norminv(1-alpha/2))/(1-a*(z0+norminv(1-alpha/2))));
		end
		bd = sort(bootdist(:,h));
		CILow = bd(round(CILowBCa*nResamples)+1);	    % lower limit of CI (following method of Wilcox 2003, see Algina et al. 2006)
		if nResamples - round((1-CIUpBCa)*nResamples) < 1
		    CIHigh = bd(1);
		else
		    CIHigh = bd(nResamples - round((1-CIUpBCa)*nResamples));  % upper limit of CI 
		end
		y(h,1) = i1;
		y(h,2) = i2;
		if j > nTraits 
		    y(h,3) = -1; 
		else 
		    y(h,3) = j;
		end
		y(h,4) = sampleValue(h);
		if j > nTraits
		    y(h,5) = 2*sampleValue(h) - mean(bootdist(:,h)); % bias correction for point estimate
		else 
		    y(h,5) = y(h,4);
		end
		y(h,6) = CILow;
		y(h,7) = CIHigh;
		y(h,8) = pValue(h); 
	    end
	end
    end
end

% -------------------------------------------------------------------------------------------------------
% function y = deltaP_single_CI(inputData, nResamples = 1000, alpha = 0.05, method = 3 useMedian = true)
% -------------------------------------------------------------------------------------------------------
% Confidence intervals for deltaP (version for single trait);
% inputData must be a matrix with trait values in the first column and
% population labels in the second column;
% nResamples is the number of bootstrap resamples;
% 1 - alpha is the confidence level;
% useMedian indicates whether deltaP is calculated based on population
% medians (default) or means;
% method indicates the method used for constructing CIs:
% method == 1: percentile method
% method == 2: BC method (bias-corrected)
% method == 3: BCa method (bias-corrected and accellerated; default)
% The output is a matrix with rows corresponding to the individual
% comparisons:
% column 1: label of population 1
% column 2: label of population 2
% column 3: point estimate for deltaP
% column 4: lower confidence limit (for the chosen method)
% column 5: upper confidence limit
% column 6: p-value from permutation test
function y = deltaP_single_CI(inputData, nResamples, alpha, method, useMedian)
    if nargin < 5
        useMedian = true;
    end
    if nargin < 4
        method = 3;
    end
    if nargin < 3
        alpha = 0.05;
    end
    if nargin < 2
        nResamples = 1000;
    end
    if size(inputData, 2) ~= 2
        error('Input data must be a n x 2 matrix.\n');
    end
    temp = deltaP_single(inputData, useMedian);
    sampleValue = temp(:,3)'; % deltaP's from original data
    temp = deltaP_single_permutation(inputData, nResamples, useMedian);
    pValue = temp(:,4)'; % p-values from permutation test
    data = inputData(:,1);
    popVector = inputData(:,2);
    popLabels = unique(popVector);
    nPops = length(popLabels);
    % create bootstrap distribution
    samplesizes = zeros(1,nPops);
    for i = 1:nPops
        sample{i} = data(popVector==popLabels(i));
        samplesizes(i) = length(sample{i});
    end
    bootdist = zeros(nResamples, nPops*(nPops-1)/2);
    for j = 1:nResamples
	resample = [];
	for i = 1:nPops
	    boot = ceil(samplesizes(i) * rand(1, samplesizes(i)));
	    clear s;
	    s = sample{i}(boot);
	    s(:,2) = popLabels(i);
	    resample = [resample; s];
	end
	temp = deltaP_single(resample, useMedian);
	bootdist(j, :) = temp(:,3)';
    end
    % create jackknife distribution (BCa method only)
    if method == 3
	jackdist = zeros(length(data), nPops*(nPops-1)/2);
	for i = 1:length(data)
	    jSample = inputData;
	    jSample(i,:) = [];
	    temp = deltaP_single(jSample, useMedian);
	    jackdist(i,:) = temp(:,3)';
	end
    end
    % calculate CIs for each comparison
    y = zeros(nPops*(nPops-1)/2, 7);
    h = 0;
    for i = 1:nPops
	for j = (i+1):nPops
	    h = h + 1;
	    % acceleration (BCa method only)
	    if method == 3
		meanJack = mean(jackdist(:,h));
		if (var(jackdist(:,h)) == 0)
		    a = 0;
		else
		    a = (sum((meanJack-jackdist(:,h)).^3))/(6*sum((meanJack-jackdist(:,h)).^2).^(3/2)); 
		end
	    else
		a = 0;
	    end
	    % bias correction (BC and BCa methods only)
	    if method >= 2
		z0 = norminv( ( sum(bootdist(:,h) < sampleValue(h)) + sum(bootdist(:,h) == sampleValue(h))/2. ) / nResamples); % bias correction
	    else
		z0 = 0;
	    end
	    if z0 == Inf 
		CILowBCa = 1;
		CIUpBCa = 1;
	    elseif z0 == -Inf
		CILowBCa = 0;
		CIUpBCa = 0;
	    else
		CILowBCa = normcdf(z0 + (z0+norminv(alpha/2))/(1-a*(z0+norminv(alpha/2))));
		CIUpBCa  = normcdf(z0 + (z0+norminv(1-alpha/2))/(1-a*(z0+norminv(1-alpha/2))));
	    end
	    bd = sort(bootdist(:,h));
	    CILow = bd(round(CILowBCa*nResamples)+1);	    % lower limit of CI (following method of Wilcox 2003, see Algina et al. 2006)
	    CIHigh = bd(nResamples - round((1-CIUpBCa)*nResamples));  % upper limit of CI 
	    y(h,1) = i;
	    y(h,2) = j;
	    y(h,3) = sampleValue(h);
	    y(h,4) = CILow;
	    y(h,5) = CIHigh;
	    y(h,6) = pValue(h);
	end
    end
end

% ------------------------------------------------------------
% function [y, P] = deltaP_single(inputData, useMedian = true)
% ------------------------------------------------------------
% DeltaP distance metric for a single trait
% inputData must be a matrix with trait values in column 1 and
% population labels in column 2; if useMedian = true, deltaP is based
% on percentiles for the population medians (default), otherwise for the
% population means
% Output: y is a nComp x 3 matrix, where nComp is the number of pairwise
% comparisons 
% Column 1: Label of population 1 in the comparison
% Column 2: Label of population 2 in the comparison
% Column 3: DeltaP for the given comparison
% P is the vector of percentiles of the joint cdf corresponding to the
% medians of the individual populations 
function [y, P] = deltaP_single(inputData, useMedian)
    if nargin < 2
        useMedian = true;
    end
    data = inputData(:,1);
    [data, index] = sort(data);
    popVector = inputData(index,2); % population of origin
    %
    popLabels = unique(popVector'); % set of population indices
    nPops = length(popLabels); % number of populations
    samplesizes = zeros(1, nPops); % sample size for each population
    center = zeros(1, nPops); % central tendencies (mean or median) of populations
    P = zeros(1, nPops); % P as in DeltaP
    for i = 1:nPops
	sample{i} = data(popVector==popLabels(i));
	samplesizes(i) = length(sample{i});
	if useMedian
	    center(i) = median(sample{i});
	else
	    center(i) = mean(sample{i});
	end
    end
    for i = 1:nPops
	for j = 1:nPops
	    P(i) = P(i) + 100 * sum(sample{j} <= center(i)) / samplesizes(j) / nPops;
	    % % In the following alternatives, data points equal to
	    % % the mean or median are only counted with a weight of 50%
	    % P(i) += 100 * (sum(sample{j} < center(i)) + sum(sample{j} == center(i))/2) / samplesizes(j) / nPops;
	end
    end
    y = zeros(nPops*(nPops-1)/2, 3);
    h = 0; % counter for ouput
    for i = 1:nPops
	for j = (i+1):nPops
	    h = h + 1;
	    y(h,1) = i;
	    y(h,2) = j;
	    y(h,3) = P(i) - P(j);
	end
    end
end

% ----------------------------------------------------------------------
% function y = deltaP_multi(inputData, useMedian = true)
% ----------------------------------------------------------------------
% DeltaP distance metric for multiple traits
% input data format: values in column 1, population labels in column 2, trait
% labels in column 3;
% If useMedian = true, deltaP is based on percentiles for the
% population medians (default), otherwise for the population means;
% Output: nComp x 4 matrix, where nComp is the number of pairwise
% comparisons 
% Column 1: Label of population 1 in the comparison
% Column 2: Label of population 2 in the comparison
% Column 3: Trait label; -1 indicates the multivariate comparison (Euclidean distance of percentiles)
% Column 4: DeltaP for the given comparison
% The multivariate comparisons appear always at the end of the output
% (thus, they can be identified even if -1 is also used as a trait
% label)
function [y,P] = deltaP_multi(inputData, useMedian)
    if nargin < 2
        useMedian = true;
    end
    if (size(inputData, 2) ~= 3)
        error('Input data must be a n x 3 matrix.\n');
    end
	traitVector = inputData(:,3);
	traitLabels = unique(traitVector);
	nTraits = length(traitLabels);
    y = [];
    P = [];
    for i = 1:nTraits
        x = inputData(traitVector == traitLabels(i), 1:2);
        [dps, ps] = deltaP_single(x, useMedian);
        dps(:,4) = dps(:,3);
        dps(:,3) = traitLabels(i) * ones(size(dps,1),1);
        y = [y; dps];
        %disp(size(P))
        %disp(size(ps))
        if i > 1  && ( size(ps,2) ~= size(P,2) )
            error('Warning:  it looks like data on at least one trait may be missing for at least one population');
        end
        P = [P; ps];
    end
    % Euclidean distances
    nPops = length(ps);
    ye = zeros(nPops*(nPops-1)/2,4);
    h = 0; 
    for i = 1:nPops
        for j = (i+1):nPops
            h = h + 1;
            ye(h,1) = i;
            ye(h,2) = j;
            ye(h,3) = -1;
            ye(h,4) = norm(P(:,i)-P(:,j));
        end
    end
    y = [y; ye];
end

% --------------------------------------------------------------------------------------
% function y = deltaP_single_permutation(inputData, nResamples = 1000, useMedian = true)
% --------------------------------------------------------------------------------------
% Two-sided permutation test of the null hypothesis deltaP = 0, for the case of
% a single trait.
% inputData must be a matrix with data values in the first column and
% population labels in the second column
% nResamples is the number of permutations
% Output: y is a nComp x 3 matrix, where nComp is the number of pairwise comparisons 
% Column 1: Label of population 1 in the comparison
% Column 2: Label of population 2 in the comparison
% Column 3: deltaP
% column 4: p-Value 
function y = deltaP_single_permutation(inputData, nResamples, useMedian)
    if nargin < 3
        useMedian = true;
    end
    if nargin < 2
        nResamples = 1000;
    end
    if (size(inputData, 2) ~= 2)
	error('Input data must be a n x 2 matrix.\n');
    end
    temp = deltaP_single(inputData, useMedian);
    sampleValue = temp(:,3)'; % deltaP's from original data
    data = inputData(:,1);
    popVector = inputData(:,2);
    popLabels = unique(popVector);
    nPops = length(popLabels);
    nComp = nPops*(nPops-1)/2;
    permdist = zeros(nResamples, nComp);
    % create and evaluate permutated data
    for k = 1:nResamples
	r = rand(size(data, 1), 1);
	[~, permut] = sort(r);
	permData = [data(permut, :), popVector];
	temp = deltaP_single(permData, useMedian);
	permdist(k, :) = temp(:,3)';
    end
    % analyze results
    y = zeros(nComp, 4);
    h = 0;
    for i1 = 1:nPops
	for i2 = (i1+1):nPops
	    h = h + 1;
	    y(h,1) = i1;
	    y(h,2) = i2;
	    y(h,3) = sampleValue(h);
	    y(h,4) = (sum(permdist(:,h) <= -abs(sampleValue(h))) + sum(permdist(:,h) >= abs(sampleValue(h)))) / nResamples;
	end
    end
    y(:,4) = min(y(:,4),1); % this should only be relevant if sampleValue == 0
end

% -------------------------------------------------------------------------------------
% function y = deltaP_multi_permutation(inputData, nResamples = 1000, useMedian = true)
% -------------------------------------------------------------------------------------
% Permutation test of the null hypothesis deltaP = 0, for data with
% multiple traits.
% input data format: values in column 1, population labels in column 2,
% trait labels in column 3
% nResamples is the number of permutations
% Output: y is a nComp x 3 matrix, where nComp is the number of pairwise comparisons 
% Column 1: Label of population 1 in the comparison
% Column 2: Label of population 2 in the comparison 
% Column 3: Trait label (-1 for multivariate deltaP, at the end of the ouptut)
% Column 4: deltaP
% column 5: p-Value 
function y = deltaP_multi_permutation(inputData, nResamples, useMedian)
    if nargin < 3
        useMedian = true;
    end
    if nargin < 2
        nResamples = 1000;
    end
    if (size(inputData, 2) ~= 3)
        error('Input data must be a n x 3 matrix.\n');
    end
    [temp,~] = deltaP_multi(inputData, useMedian);
    sampleValue = temp(:,4)'; % deltaP's from original data
	data = inputData(:,1);
	popVector = inputData(:,2);
	traitVector = inputData(:,3);
	traitLabels = unique(traitVector);
	nTraits = length(traitLabels);
    popLabels = unique(popVector);
    nPops = length(popLabels);
    nComp = (nTraits+1)*nPops*(nPops-1)/2;
    permdist = zeros(nResamples, nComp);
    % create and evaluate permutated data
    for k = 1:nResamples
	    permData = [];
	    for i = 1:nTraits
            s = data(traitVector==traitLabels(i));
            r = rand(length(s),1);
            [~,permut] = sort(r);
            s = s(permut);
            s(:,2) = popVector(traitVector==traitLabels(i));
            s(:,3) = i;
            permData = [permData; s];
        end
        [temp,~] = deltaP_multi(permData, useMedian);
        permdist(k, :) = temp(:,4)';
    end
    % analyze results
    y = zeros(nComp, 5);
    h = 0;
    for j = 1:nTraits+1
	for i1 = 1:nPops
	    for i2 = (i1+1):nPops
		h = h + 1;
		y(h,1) = i1;
		y(h,2) = i2;
		if j <= nTraits 
		    y(h,3) = j; 
		    y(h,5) = (sum(permdist(:,h) <= -abs(sampleValue(h))) + sum(permdist(:,h) >= abs(sampleValue(h)))) / nResamples;
		else 
		    y(h,3) = -1;
		    y(h,5) = sum(permdist(:,h) >= abs(sampleValue(h))) / nResamples;
		end
		y(h,4) = sampleValue(h);
	    end
	end
    end
    y(:,5) = min(y(:,5),1); % this is only relevant if sampleValue == 0
end



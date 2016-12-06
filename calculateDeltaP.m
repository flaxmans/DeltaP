% script to import data from a plain text file and run calculations of
% deltaP
%
% updated 10 Feb. 2012

% HOW TO USE THE SOURCE CODE TO CALCULATE DELTA-P and CONFIDENCE INTERVALS:  
% (1) this file plus the function files called below
% (i.e., "deltaPinput.m", "deltaP_multi_CI.m", etc.)
% must all be together in one directory, which must
% be your current working directory (or part of your path) when you execute this
% (2) either type "calculateDeltaP" at the MATLAB command line, or simply
% press MATLAB's "Run" button from the Editor window
% (3) Respond to user input
% (4) Report bugs/failures to Sam Flaxman (samuel.flaxman@colorado.edu)

% if the user chooses, the following values can be changed; if not, all
% will run with defaults

nResamples = 1000;  % number of resamples in bootstrapping procedure
alpha = 0.05;  % alpha level for construction of confidence intervals
useMedian = true;  % use median of each observation set for percentile calculation (alternative is mean if this is set to false)
method = [3,1]; % bootstrapping procedure; see comments in deltaP_multi_CI.m for explanation


% import the data
[~,~,~,~,~,~,inputData,traitNames,populationNames,fileName] = deltaPinput();
% see comments at the top of the deltaPinput.m file if the import of data 
% is not working properly


% calculate deltaP values and associated confidence intervals
[deltaPs, rawPercentiles] = deltaP_multi_CI(inputData, nResamples, alpha, method, useMedian);
README_deltaPs={'column 1: index of first population, i.e., populationNames(deltaPs(:,1))';...
    'column 2: index of second population, i.e., populationNames(deltaPs(:,2))';...
    'column 3: trait index; given as "-1" for the euclidean distances';...
    'column 4: point estimate for deltaP';...
    'column 5: bias-corrected point estimate; only applied for the multivariate eucl. dist.';...
    'column 6: lower confidence limit (for the chosen method)';...
    'column 7: upper confidence limit';...
    'column 8: p-value from permutation test'};


% calculate Hedges' G statistics for all pairwise comparisons and all
% traits
Hedges_G = calculateHedgesG(inputData);
rows = size(Hedges_G,1);
Hedges_G = [deltaPs(1:rows,1:3) Hedges_G];
README_Hedges_G = {'column 1: index of first population, i.e., populationNames(Hedges_G(:,1))';...
    'column 2: index of second population, i.e., populationNames(Hedges_G(:,2))';...
    'column 3: trait index';...
    'column 4: Hedges G statistic'};
clear('rows');


% write the output to a file
outputDeltaPs(deltaPs,Hedges_G(:,4),rawPercentiles,inputData,traitNames,populationNames,fileName);
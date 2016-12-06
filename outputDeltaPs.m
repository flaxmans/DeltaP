function outputDeltaPs(dps,G,rps,data,traitNames,popNames,fileName)
% "dps" is the n x 6 matrix of delta-P values returned from the function
% dps = deltaP_multi_CI() called above.  This matrix is all numbers
% col 1: population 1 (as integer)
% col 2: population 2 (as integer)
% col 3: trait (as integer)
% col 4: deltaP (point estimate)
% col 5: the bias-corrected deltaP mentioned in Appendix 4 of Safran et al. 2012
% col 6: lower boundery of CI
% col 7: upper boundery of CI
% col 8: p-value from permutation test


% 1st column of dps = pop. 1 (as integer)
% 2nd column = pop. 2 (as integer)
% 3rd column = trait (as integer)
% 4th column = point estimate of delta-P
% 5th column = lower bound of CI
% 6th column = upper bound of CI
% rps = raw percentiles: ntraits x npops matrix

disp(' ');
disp('Results on delta-P will be output to the following file:  ');
outName = [fileName '_deltaPs.csv'];
disp(['    ' outName]);
disp('Note: continuing will overwrite the file of the same name (if it exists).');
answer = input('Please press return to continue, or enter a different file name:  ','s');
disp(' ');
if strcmp(answer,'')
    fpt = fopen(outName,'w');
    disp(['   *** Writing output to "' outName '" ***']);
else
    fpt = fopen(answer,'w');
    disp(['   *** Writing output to "' answer '" ***']);
end
disp(' ');

%traits one at a time
npops = length(popNames);
meaningfulComps = npops * (npops - 1) / 2;
si = 1;
ei = meaningfulComps;
for i = 1:length(traitNames)
    fprintf(fpt,'Results for trait %i: %s\n',i,traitNames{i});
    
    descriptiveStats = traitSummaries(data,i);
    outputDS(descriptiveStats,fpt,popNames,rps(i,:));
    
    fprintf(fpt,'Popn. 1,Popn. 2,delta-P,lower CI limit,upper CI limit,p-value from permutation test,Hedges G\n');
    [~,n] = sort(abs(dps(si:ei,4)),'descend');
    for j = 1:meaningfulComps
        wr = si - 1 + n(j); % working row based upon sort
        fprintf(fpt,'%s,%s,%f,%f,%f,%f,%f\n', popNames{dps(wr,1)},popNames{dps(wr,2)},dps(wr,4),dps(wr,6),dps(wr,7),dps(wr,8),G(wr));
    end
    si = ei + 1;
    ei = ei + meaningfulComps;
    fprintf(fpt,'\n');
end

%summary of population distances
[~,n] = sort(dps(si:ei,4),'descend');
fprintf(fpt,'Summary of population distances from each other based upon all traits\n');
fprintf(fpt,'Popn. 1,Popn. 2,Euclidean Distance,bias-corrected E.D.,lower CI limit,upper CI limit,p-value from permutation test\n');
for j = 1:meaningfulComps
    wr = si - 1 + n(j); % working row based upon sort
    fprintf(fpt,'%s,%s,%f,%f,%f,%f,%f\n', popNames{dps(wr,1)},popNames{dps(wr,2)},dps(wr,4),dps(wr,5),dps(wr,6),dps(wr,7),dps(wr,8));
end

% raw percentiles
fprintf(fpt,'\nPercentiles for all traits and populations (concatenation of "percentile" columns above)\n');
for i = 1:npops
    fprintf(fpt,',%s',popNames{i});
end
fprintf(fpt,'\n');
for i = 1:length(traitNames)
    fprintf(fpt,'%s',traitNames{i});
    for j = 1:npops
        fprintf(fpt,',%f',rps(i,j));
    end
    fprintf(fpt,'\n');
end
        





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate descriptive statistics
function traitStats = traitSummaries(input,ti)
npops = length(unique(input(:,2)));
traitStats = zeros(npops,7);
% 7 columns of traitStats will be, in order: sample size, min, max, mean,
% median, sd, cv
wt = find(input(:,3)==ti); % rows containing trait of interest
for i=1:npops
    wp = find(input(:,2) == i);
    wis = intersect(wt,wp);
    vals = input(wis,1);
    traitStats(i,1) = length(vals);
    traitStats(i,2) = min(vals);
    traitStats(i,3) = max(vals);
    traitStats(i,4) = mean(vals);
    traitStats(i,5) = median(vals);
    traitStats(i,6) = std(vals);
    traitStats(i,7) = traitStats(i,6) / traitStats(i,4);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to write descriptive statistics to output file
function outputDS(DS,fpt,popNames,rps)
fprintf(fpt,'population,sample size,min,max,mean,median,std. dev.,c.v.,percentile\n');
rows = size(DS,1);
for i = 1:rows
    fprintf(fpt,'%s',popNames{i});
    for j = 1:7
        fprintf(fpt,',%f',DS(i,j));
    end
    fprintf(fpt,',%f\n',rps(i));
end
fprintf(fpt,'\n');
    



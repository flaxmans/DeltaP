function [r,traitN,traitT,populationN,populationT,values,ValPopTraitMatrix,traitNames,popNames,fileName] = deltaPinput()

% deltaPinput.m
% last revised 9 February 2012
% Copyright 2011-2012, Samuel Melvin Flaxman and Michael Kopp, all rights reserved
% please report bugs or errors to samuel.flaxman@colorado.edu


% function to import data on multiple traits from multiple populations

% The function assumes that the data are in a format that Matlab
% understands, such as a comma separated values format

% Matlab's importdata() function is fairly robust to missing cells, but
% this will be most reliable if nearly all cells are filled and if a column
% of data is either ALL text or ALL numbers (with the exception of the
% header row).  Mixed columns or extensive gaps in the data may prevent a
% proper import.  The function as written below does a little self
% reporting on what it has imported, and this should be checked against
% what one sees in the raw data file

% The function as written assumes implicitly that the data file is in the
% same directory as the function OR that you will enter the "full name"  
% as ${PATH}/filename , including the extension in the filename









%%%%%%%%%%%    USER INPUT %%%%%%%%%%%%%%%

% much of this code is here to try to handle the vagaries of real data and
% human error/inconsistency in recording and formatting it


% first, get the data
[r,nhead,proceed,rismatrix,risstruct,fileName] = initialImportOfData();
if ~proceed
    return;
end


if risstruct % matlab made the import into a structure with text and numerical data
    % next, create vectors that store the indexes of different categories of
    % columns
    textColumns = zeros(1,size(r.textdata,2) - size(r.data,2));
    % vector of indexes of columns imported as text by MATLAB
    dataColumns = zeros(1,size(r.data,2));
    % vector of indexes of columns imported as numbers by MATLAB
    nextt = 1;
    nextd = 1;
    for i = 1:size(r.textdata,2)
        if strcmp('',r.textdata(2,i)) % number
            dataColumns(nextd) = i;
            nextd = nextd + 1;
        else % text
            textColumns(nextt) = i;
            nextt = nextt + 1;
        end
    end
elseif rismatrix  % matlab imported only numbers
    dataColumns = 1:nhead;
    textColumns = [];
else
    disp(' ');
    disp('    ERROR!  imported data not in recognized format!  Exiting.');
    return;
end

% get the population names in vectors as strings and numbers
[populationT,populationN,proceed] = parsePopulationIDs(textColumns,dataColumns,r,rismatrix,risstruct);   
if ~proceed
    disp(' ');
    disp('    ERROR in read of population names.  Please check data file.');
    disp(' ');
    return;
end

% next, deal with stacked or unstacked trait data
disp(' ');
answer = input('Are you importing data on more than one trait?  ','s');
multTraits = resolveYesNoAnswer(answer);
if multTraits
    answer = input('Are data on MULTIPLE traits stacked in a SINGLE column of data?  ','s');
    stacked = resolveYesNoAnswer(answer);
else
    stacked = 0;
end

if stacked
    [traitT,traitN,values,proceed] = parseStackedTraitData(textColumns,dataColumns,r,rismatrix,risstruct);
else
    [traitT,traitN,values,traitNames,populationT,populationN,proceed] = parseUnstackedTraitData(dataColumns,r,rismatrix,risstruct,populationT,populationN,multTraits);
end
if ~proceed
    disp(' ');
    disp('    ERROR in parsing of trait names.  Please check data file.');
    disp(' ');
    return;
end

disp(' ');
answer = input('How many populations or distinct groups were sampled?  ');
expectedN = resolveIntegerAnswer(answer);
if expectedN < 2
    disp(' ');
    disp('    ERROR!  At least two populations are required for this method.');
    disp('    exiting function');
    disp(' ');
    return;
end
[popNames,npops,emptyCells] = parseCellArray(populationT,expectedN);
if npops ~= expectedN
    disp(' ');
    disp('   ERROR!  number of populations does not equal expected number!');
    disp('   Check data file for typos and/or missing/additional populations');
    disp('   Exiting.');
    return;
elseif emptyCells
    disp(' ');
    disp(['   ERROR!  found ' num2str(emptyCells) ' empty cells in population names column!']);
    disp('   Check data file for missing names.');
    disp('   Exiting.');
    return;
end
disp(' ');
disp('Found the following population names:');  % report what was found
disp(popNames);

if stacked
    disp(' ');
    answer = input('How many traits were measured?  ');
    expectedN = resolveIntegerAnswer(answer);
    [traitNames,numTraits,emptyCells] = parseCellArray(traitT,expectedN);
    if numTraits ~= expectedN
        disp(' ');
        disp('   ERROR!  number of traits does not equal expected number!');
        disp('   Check data file for typos and/or missing/additional trait names');
        disp('   Exiting.');
        return;
    elseif emptyCells
        disp(' ');
        disp(['   ERROR!  found ' num2str(emptyCells) ' empty cells in trait names column!']);
        disp('   Check data file for missing names.');
        disp('   Exiting.');
        return;
    end
end

disp(' ');
disp('Found the following trait names:');  % report what was found
disp(traitNames);


ValPopTraitMatrix = [values populationN traitN];



%%%%%%%%%%% END OF USER INPUT FOR DATA IMPORT %%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  HERE AND BELOW:  Accessory functions used for I/O operations %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,nhead,proceed,rismatrix,risstruct,targetFile] = initialImportOfData()
disp(' ');
disp('Enter full name (including file extension, if applicable)');
targetFile = input('of data file to import:  ','s');

disp(' ');
disp('    *** Looking for File ***');

r = importdata(targetFile);

disp(' ');
disp('    *** Importing Data ***');
disp(' ');

if iscell(r)
    disp(' ');
    disp('It appears that the imported data were interpreted as strings');
    disp('and not as numbers.  Please reformat the data file so that text');
    disp('identifiers are all in the left-most columns and numerical');
    disp('values are all in the right-most columns.');
    disp(' ');
    input('Press return to display imported data and exit function.');
    disp(' ');
    disp('imported data = ');
    disp(r);
    disp(' ');
    rismatrix = 0;
    risstruct = 0;
    proceed = 0;
    return;
end
%format long;

if isstruct(r)
    risstruct = 1;
    rismatrix = 0;
    [datarows nhead] = size(r.textdata);
    if datarows == 1 % r.textdata is ONLY the column headers
        dumCell = cell(size(r.data));
        for i = 1:numel(dumCell)
            dumCell{i} = '';
        end
        r.textdata = [r.textdata; dumCell];
        datarows = size(r.textdata,1);
    end
    datarows = datarows - 1; % number of actual rows of data (not including header)
    if ( datarows ~= size(r.data,1) ) % checking for good import of file meeting assumptions
        disp('Error in import!');
        if datarows == size(r.data,1) - 1
            disp('Data file appears to have no text header row.')
            disp('Exiting function');
            disp(' ');
        elseif datarows < size(r.data,1)
            disp('Check for blank rows at bottom of data sheet or text characters mixed with numbers.')
            disp(' ');
            disp('Exiting function.  Try running this function again after revising and saving data sheet.');
            disp(' ');
        end
        return;
    end

    disp(' ');
    disp('    *** Import Completed ***');
    disp(' ');


    disp(' ');
    disp(['Data file imported with ' num2str(datarows) ' rows of data, ']);
    disp([num2str(size(r.textdata,2) - size(r.data,2)) ' columns of text data,']);
    disp([num2str(size(r.data,2)) ' columns of numerical data,']);
    disp(['and the following ' num2str(nhead) ' column headers:  ']);
    disp(' ');
    disp(r.textdata(1,:));
    disp('Does this look correct?');
    disp('Enter "yes" to proceed; enter "no" to display what was imported');
    answer = input('and to exit the function:  ','s');
    proceed = resolveYesNoAnswer(answer);
    if ~proceed
        disp('text data = ');
        disp(r.textdata);
        disp('numerical data = ');
        disp(r.data);
    end
elseif isnumeric(r) % the import was just columns of numbers
    [datarows,nhead] = size(r);
    rismatrix = 1;
    risstruct = 0;
    proceed = 1;
    disp(' ');
    disp('    *** Import Completed ***');
    disp(' ');
    disp(['Data file imported with ' num2str(datarows) ' rows of data, ']);
    disp(['and ' num2str(nhead) ' columns of numerical data, and no column headers.']);
    disp(' ');
end


    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [populationT,populationN,proceed] = parsePopulationIDs(textColumns,dataColumns,r,rismatrix,risstruct)

disp('Which column identifies the population from which an observation was taken?');
answer = input('Enter the column number as an integer:  ');
popCol = resolveIntegerAnswer(answer);

if risstruct
    
    rowsOfData = size(r.data,1);
    if any(dataColumns == popCol)
        popColWithinData = find(dataColumns == popCol);
        populationN = r.data(:,popColWithinData);
        populationT = cell(rowsOfData,1);
        for i = 1:rowsOfData
            populationT{i} = num2str(r.data(i,popColWithinData));
        end
    elseif any(textColumns == popCol)
        popColWithinText = find(textColumns == popCol);
        populationT = r.textdata(2:(rowsOfData+1),popColWithinText);
        populationN = zeros(rowsOfData,1);
        popNames = cell(rowsOfData,1);
        found = 1;
        popNames{1} = populationT{1};
        populationN(1) = 1;
        for i = 2:rowsOfData
            if any(strcmp(populationT{i},popNames))
                populationN(i) = find(strcmp(populationT{i},popNames));
            else
                found = found + 1;
                populationN(i) = found;
                popNames{found} = populationT{i};
            end
        end     
    else
        disp(' ');
        disp('     ****** ERROR!  Selected column appears to be invalid! *******');
        disp(' ');
        proceed = 0;
        return;
    end 
    
elseif rismatrix
    
    rowsOfData = size(r,1);
    populationN = r(:,popCol);
    populationT = cell(size(populationN));
    for i = 1:rowsOfData
        populationT{i} = num2str(populationN(i));
    end
    
end
proceed = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [traitT,traitN,values,proceed] = parseStackedTraitData(textColumns,dataColumns,r,rismatrix,risstruct)

disp('Which column identifies the trait measured for a given observation value?');
answer = input('Enter the column number as an integer:  ');
traitsCol = resolveIntegerAnswer(answer);
if risstruct

    rowsOfData = size(r.data,1);
    if any(dataColumns == traitsCol)
        traitsColWithinData = find(dataColumns == traitsCol);
        traitN = r.data(:,traitsColWithinData);
        traitT = cell(rowsOfData,1);
        for i = 1:rowsOfData
            traitT{i} = num2str(r.data(i,traitsColWithinData));
        end
    elseif any(textColumns == traitsCol)
        traitsColWithinText = find(textColumns == traitsCol);
        traitT = r.textdata(2:(rowsOfData+1),traitsColWithinText);
        traitN = zeros(rowsOfData,1);
        traitNames = cell(rowsOfData,1);
        found = 1;
        traitNames{1} = traitT{1};
        traitN(1) = 1;
        for i = 2:rowsOfData
            if any(strcmp(traitT{i},traitNames))
                traitN(i) = find(strcmp(traitT{i},traitNames));
            else
                found = found + 1;
                traitN(i) = found;
                traitNames{found} = traitT{i};
            end
        end     
    else
        disp(' ');
        disp('     ****** ERROR!  Selected column appears to be invalid! *******');
        disp(' ');
        proceed = 0;
        return;
    end 
    
elseif rismatrix
    
    rowsOfData = size(r,1);
    traitN = r(:,traitsCol);
    traitT = cell(size(traitN));
    for i = 1:rowsOfData
        traitT{i} = num2str(traitN(i));
    end
    
end




disp('Which column identifies the measured value of a trait?');
answer = input('Enter the column number as an integer:  ');
valuesCol = resolveIntegerAnswer(answer);
if risstruct
    if any(dataColumns == valuesCol)
        valuesColWithinData = find(dataColumns == valuesCol);
        values = r.data(:,valuesColWithinData);
    else
        disp(' ');
        disp('     ****** ERROR!  Selected column appears to be invalid! *******');
        disp(' ');
        proceed = 0;
        return;
    end
elseif rismatrix
    values = r(:,valuesCol);
end

    
proceed = 1;

% disp(populationT);
% disp(populationN);
% disp(traitT);
% disp(traitN);
% disp(values);
% return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [traitT,traitN,values,traitNames,populationT,populationN,proceed] = parseUnstackedTraitData(dataColumns,r,rismatrix,risstruct,populationT,populationN,multTraits)

if multTraits 
    disp('Which columns give measurements of traits?');
    tcs = input('Enter column numbers as integers separated by commas:  ','s');
else
    disp('Which column gives measurements of the trait?');
    tcs = input('Enter the column number as an integer:  ','s');
end
eval(['traitCols = [' tcs '];']);
numTraits = length(traitCols);
traitNames = cell(1,numTraits);

if risstruct

    rowsOfData = size(r.data,1);
    unstackedTraits = zeros(rowsOfData,numTraits);
    traitN = zeros((numTraits*rowsOfData),1);
    traitT = cell(size(traitN));
    values = traitN;
    newPT = traitT;
    newPN = traitN;
    
    for i = 1:numTraits
        dum = find(dataColumns == traitCols(i));
        wc = r.data(:,dum); % working column of data
        unstackedTraits(:,i) = wc;
        traitName = r.textdata{1,traitCols(i)};
        traitNames{i} = traitName;
        si = 1 + ((i-1)*rowsOfData);
        ei = i * rowsOfData;
        values(si:ei) = wc;
        newPT(si:ei) = populationT;
        newPN(si:ei) = populationN;
        for j=si:ei
            traitT{j} = traitName;
            traitN(j) = i;
        end
    end
    
elseif rismatrix
    madeUpNames = 1:numTraits;
    traitNames = cell(size(madeUpNames));
    for i = 1:length(madeUpNames)
        traitNames{i} = num2str(madeUpNames(i));
    end
    rowsOfData = size(r,1);
    unstackedTraits = zeros(rowsOfData,numTraits);
    traitN = zeros((numTraits*rowsOfData),1);
    traitT = cell(size(traitN));
    values = traitN;
    
    for i = 1:numTraits
        wc = r(traitCols(i)); % working column of data
        traitName = traitNames{i};
        si = 1 + ((i-1)*rowsOfData);
        ei = i * rowsOfData;
        values(si:ei) = wc;
        unstackedTraits(:,i) = wc;  % this isn't used in this version, but is still calculated in case a need arises
        newPT(si:ei) = populationT;
        newPN(si:ei) = populationN;
        for j=si:ei
            traitT{j} = traitName;
            traitN(j) = i;
        end
    end
    
end
    
if any(isnan(values))
    removals = find(isnan(values));
    traitT(removals) = [];
    traitN(removals) = [];
    values(removals) = [];
    newPT(removals) = [];
    newPN(removals) = [];
end

populationT = newPT;
populationN = newPN;

proceed = 1;
        
        




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [names,found,emptyCells] = parseCellArray(array,expectedN)

emptyCells = 0;
found = 1;
names = cell(1,expectedN);
if strcmp('',array{1})
    disp('Warning!  Expected name but found empty cell!');
    emptyCells = 1;
else
    names{1} = array{1};
    for i = 2:length(array)
        if strcmp('',array{i})
            disp('Warning!  Expected name but found empty cell!');
            emptyCells = emptyCells + 1;
        elseif ~any(strcmp(array{i},names))
            found = found + 1;
            if found > expectedN
                disp('Warning!  Found more names than expected!  Check data file for typos or additional names!');
                disp(['extra name found = ' array{i}]);
                %input('Press return to continue or CNTRL-C to exit  ');
            else
                names{found} = array{i};
            end
        end
    end
end
if found < expectedN
    disp('Warning!  Found fewer names than expected!');
    %disp('Check data file for missing names');
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = resolveYesNoAnswer(answer)
notUnderstood = 1;
yesses = {'Y','y','YES','YEs','YeS','Yes','yES','yEs','yeS','yes'};
nos = {'N','n','NO','no','No','nO'};
disp(' ');
while notUnderstood
    if any(strcmp(answer,yesses))
        bool = 1;
        notUnderstood = 0;
    elseif any(strcmp(answer,nos))
        bool = 0;
        notUnderstood = 0;
    else 
        disp('Response not understood.');
        answer = input('Please enter yes or no:  ','s');
    end
end
disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intValue = resolveIntegerAnswer(answer)
if isnumeric(answer)
    intValue = floor(answer);
else
    isNotInteger = 1;
    while isNotInteger
        intValue = input('Please enter an INTEGER:  ');
        if isnumeric(intValue)
            isNotInteger = 0;
            intValue = floor(intValue);
        end
    end
end



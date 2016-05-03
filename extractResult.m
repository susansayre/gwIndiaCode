function resultList = extractResult(resultsArray,resultString)

cases = numel(resultsArray);
resultList = [];
for ii=1:cases
    eval(['resultList = [ resultList resultsArray{ii}.' resultString '];'])
end
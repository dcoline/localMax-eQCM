%script for ttest between avg (to be PCA) of response/non-response per cluster > discard p-value <.1

%stratified LC patient samples by column number

noResponse= keepSMALLData(:,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103]);
response= keepSMALLData(:,[104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181]);

%avg expression of strata collapsed by sample
for i=1:1:length(sizeMergedCluster)
clusterData= response(mergedCluster{i}(2:end), :);
meanClusterData = mean(clusterData);
clusterData2= noResponse(mergedCluster{i}(2:end), :);
NoMeanClusterData = mean(clusterData2);

%ttest among strata avg/PCA
[h,p] = ttest2 ((meanClusterData),(NoMeanClusterData),'Alpha',0.05);

%display to array
%intraclusterPvalue=(p, strcat('cluster' {i} 'p-value', num2str(i))p-value{i+1}));
pv(i)=[p]
end
%keep clusters pv<.05
for i=1:1:length(sizeMergedCluster)
if pv(i)<.05
      intraClusterDiff = [keepGeneSym(mergedCluster{i})]
end


end
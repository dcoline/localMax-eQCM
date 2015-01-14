%script for ttest between avg (to be PCA) of response/non-response per cluster > discard p-value <.1

%stratified LC patient samples by column number

noResponse= keepSMALLData(:,[13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
response= keepSMALLData(:,[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22]);


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
if pv(i)<.1
      intraClusterDiff = [keepGeneSym(mergedCluster{i})]
end


end
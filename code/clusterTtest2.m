%script for ttest between avg (to be PCA) of response/non-response per cluster > discard p-value <.1

%stratified LC patient samples by column number

response= keepSMALLData(:,[16,91,44,73,22,71,74,72,13,58,85,102,82,14,33,3,18,25]);
noResponse= keepSMALLData(:,[60,39,81,46,12,17,32,34,19,53,95]);
% 
% 
% %avg expression of strata collapsed by sample
% for i=1:1:length(sizeMergedCluster)
% clusterData= response(mergedCluster{i}(2:end), :);
% meanClusterData = mean(clusterData);
% clusterData2= noResponse(mergedCluster{i}(2:end), :);
% NoMeanClusterData = mean(clusterData2);
% 
% %ttest among strata avg/PCA
% [h,p] = ttest2 ((meanClusterData),(NoMeanClusterData),'Alpha',0.05);
% 
% %display to array
% %intraclusterPvalue=(p, strcat('cluster' {i} 'p-value', num2str(i))p-value{i+1}));
% pv(i)=[p]
% end
% %keep clusters pv<.05
% for i=1:1:length(sizeMergedCluster)
% if pv(i)<.05
%       intraClusterDiff = [keepGeneSym(mergedCluster{i})]
% end
% 
% 
% end
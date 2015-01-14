
noResponse= dataMatrix(:,[1,28,29,56,48]);
response= dataMatrix(:,[63,8,26,52,3,57]);


%avg expression of strata collapsed by sample
for i=1:1:length(sizeMergedCluster)
temp{i}= dataMatrix(mergedCluster{i}, :);
%temp{i}=(temp{i})'
[u,s,v] = svd(temp{i});
V{i}=[v(:,1)]
end


% 
fid = fopen(strcat(filepath, 'anovapairs_Fil.txt'), 'w');
for i = 1 : length(mergedCluster)
    fprintf(fid, '%s\t', strcat('cluster ', num2str(i),'eigengene'));
    geneList = V(1,i);
    for j = 1 : length(mergedCluster{i})
        fprintf(fid, '%s\t', V{1,i});
    end;
    fprintf(fid, '\n');
end;
fclose(fid);


%ttest among strata avg/PCA
% for i=1:1:length(sizeMergedCluster)
% W{i} = V{1,i}{1,28,29,56,48};
% X{i} = V{1,i}{63,8,26,52,3,57};
% % %display to array
% % %intraclusterPvalue=(p, strcat('cluster' {i} 'p-value', num2str(i))p-value{i+1}));
%  end
% %pv2(i)=1-pv(i)
% 
% %keep clusters pv<.05
% for i=1:1:length(sizeMergedCluster)
% if p{i}<.1
%       intraClusterDiff{i} = [geneSym(mergedCluster{i})]
% end
% end
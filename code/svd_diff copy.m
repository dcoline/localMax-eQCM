%script for ttest between avg (to be PCA) of response/non-response per cluster > discard p-value <.1
%stratified LC patient samples by column number
%noResponse= keepSMALLData(:,[13,14,25,16,74,76,40,44,64,68]);
%response= keepSMALLData(:,[8,12,15,16,17,19,21,22,28,32,37,40,42,43,53,54,55,57,59,63,73,77,81,82,95]);
%response= keepSMALLData(:,[22,37,17,6,72,94,73,82,8,18]);
% 12,15,19,28,32,42,43,57,77,81,95
%%%noResponse= keepSMALLData(:,[28,53,47,44,41,99,62,13,14,72,74,76,25]);
%%%response= keepSMALLData(:,[21,63,55,54,8,82,6,40,16,73,17,37,22]);

noResponse= dataMatrix(:,[1,28,29,56,48]);
response= dataMatrix(:,[63,8,26,52,3,57]);

%stratified LC patient samples by column number
% MT1= keepSMALLData([582,1551,1222,4936,2428,3180],[3,13,14,18,25,44,64,68,72,75,76,85,94,100,102]);
% MT2= keepSMALLData([582,1551,1222,4936,2428,3180],[8,12,15,16,17,19,21,22,28,32,37,40,42,43,53,54,55,57,59,63,73,77,81,82,95]);
% HIST= keepSMALLData([798,2128,5862,5940,7857,9989],[3,13,14,18,25,44,64,68,72,75,76,85,94,100,102,8,12,15,16,17,19,21,22,28,32,37,40,42,43,53,54,55,57,59,63,73,77,81,82,95]);
% IFIT= keepSMALLData([6756,9588,9075,9323,9225],[3,13,14,18,25,44,64,68,72,75,76,85,94,100,102,8,12,15,16,17,19,21,22,28,32,37,40,42,43,53,54,55,57,59,63,73,77,81,82,95]);
% LOC= keepSMALLData([9440,10514,12526],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
% ATA= keepSMALLData([1110,1111],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
% NFIA= keepSMALLData([12308,8442],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
%IFITM= data([7749,8442,8110,1054,9434,1324,8439,1994],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);

%avg expression of strata collapsed by sample
for i=1:1:length(sizeMergedCluster)
temp{i}= response(mergedCluster{i}, :);
%temp{i}=(temp{i})'
[u,s,v] = svd(temp{i});
yes{i}=[v(:,1)]

temp2{i}= noResponse(mergedCluster{i}, :);
%temp2{i}=(temp2{i})'
[nu,ns,nv] = svd(temp2{i});
no{i}=[nv(:,1)]

g{i}=yes{i},no{i};
end
fid = fopen(strcat(filepath, 'anovapairs3.txt'), 'w');
for i = 1 : length(mergedCluster)
    fprintf(fid, '%s\t', strcat('pair ', num2str(i)));
    geneList = g(1,i);
    for j = 1 : length(mergedCluster{i})
        fprintf(fid, '%s\t', g{1,i}{1,2});
    end;
    fprintf(fid, '\n');
end;
fclose(fid);

for i=1:1:length(sizeMergedCluster)
p{i} = anova1 (g{1,i}{1,1}),(g{1,i}{1,2});

%display to array
%intraclusterPvalue=(p, strcat('cluster' {i} 'p-value', num2str(i))p-value{i+1}));
end
%pv2(i)=1-pv(i)

%keep clusters pv<.05
for i=1:1:length(sizeMergedCluster)
if pv(i)<.1
      intraClusterDiff{i} = [geneSym(mergedCluster{i})]
end
end
% %%%%%%%%%%% Step 1 - Compute PCC matrix %%%%%%%%%%%%%%%%
% 
% %%%% if there is NaN or missing data in the data matrix, use below %%%%%%
% tic
% cMatrix = zefros(nGene, nGene);
% for i = 1 : nGene-1
%     if (mod(i, 200) == 0)
%         i
%         toc
%     end;
%     cMatrix(i+1:end,i) = massivePCC_withNaN(sampleData(i+1:end, :), sampleData(i,:), 10);
% end;
% cMatrix = cMatrix + cMatrix';
% 
% %%%%% if there is NO NaN or missing data in the data matrix, use below %%%
% cMatrix = massivePCC_withoutNaN(sampleData);
cMatrix(1 : size(cMatrix, 1)+1 : end) = 0;
%%%%%%%% if weigth normalization is needed, use below %%%%%%%%%%%%%%%%%%%%
cMatrix = abs(cMatrix);
D = sum(cMatrix);
D_half = 1./sqrt (D);
for i = 1 : size(cMatrix, 1)
    cMatrix(i, :) = cMatrix(i, :) * D_half(i);
end;
for i = 1 : size(cMatrix, 1)
    cMatrix(:, i) = cMatrix(:, i) * D_half(i);
end;


%%%%%%%% Step 2 - identify co-expression modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Algorithm parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = .4; t = 1; lambda = 1;
%%%%%%%% Run the algorithm
C = localMaximumQCM(abs(cMatrix), gamma, t, lambda);
%%%%%%%% Step 3 - Merge the overlapped networks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Allowed overlap ratio threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.4; 

%%%%%%%% Sort and iteratively merge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeC = zeros(1, length(C));
for i = 1 : length(C)
    sizeC(i) = length(C{i});
end;

[sortC, sortInd] = sort(sizeC, 'descend');                                  %step 4.1
C = C(sortInd);

ind = find(sortC >= 3);

mergedCluster = C(ind);
mergeOccur = 1; 
currentInd = 0;

while mergeOccur == 1
    mergeOccur = 0;
    while currentInd < length(mergedCluster)
        currentInd = currentInd + 1;                                        %step 4.4
        excludeInd = [];
        if (currentInd < length(mergedCluster))
            keepInd = 1 : currentInd;
            for j = currentInd+1 : length(mergedCluster)                    %step 4.3
                interCluster = intersect(mergedCluster{currentInd}, mergedCluster{j});
                if length(interCluster) >= beta*min(length(mergedCluster{j}), length(mergedCluster{currentInd}))%step 4.2
                    mergedCluster{currentInd} = union(mergedCluster{currentInd}, mergedCluster{j});
                    mergeOccur = 1;
                else
                    keepInd = [keepInd, j];
                end;
            end;
            mergedCluster = mergedCluster(keepInd);
            length(mergedCluster)
        end;
    end;
    sizeMergedCluster = zeros(1, length(mergedCluster));
    for i = 1 : length(mergedCluster)
        sizeMergedCluster(i) = length(mergedCluster{i});
    end;
    [sortSize, sortMergedInd] = sort(sizeMergedCluster, 'descend');
    mergedCluster = mergedCluster(sortMergedInd);
    currentInd = 0;
end;

fid = fopen(strcat(filepath, 'clusterGeneList_FILT_NORM.txt'), 'w');
for i = 1 : length(mergedCluster)
    fprintf(fid, '%s\t', strcat('cluster ', num2str(i)));
    geneList = geneSym(mergedCluster{i});
    for j = 1 : length(geneList)
        fprintf(fid, '%s\t', geneList{j});
    end;
    fprintf(fid, '\n');
end;
fclose(fid);
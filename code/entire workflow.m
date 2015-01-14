clear;
close all;

LOAD

filepath = '/Users/danielcolin/Dropbox/OSU/Spring14/Thesis KUN/TWO/code/';
datapath = '/Users/danielcolin/Dropbox/OSU/Spring14/Thesis KUN/TWO/code/';
datafile = 'ccleLCexpr2.txt';
fid = fopen(strcat(datapath, datafile));
tmpLine = fgetl(fid); 
A = textscan(tmpLine, '%s', 'delimiter', '\t'); A = A{1};
nSamples = length(A) - 2;
textstr = strcat('%s %s ', repmat('%f ', 1, nSamples));
B = textscan(fid, textstr, 'delimiter', '\t');
fclose(fid);

geneSym = B{1}; probeId = B{2};
dataMatrix = [];
for i = 1 : nSamples
    dataMatrix = [dataMatrix, B{i+2}];
end;
save(strcat(datapath, 'CCL_LC_Parsed.txt'), 'geneSym', 'dataMatrix', 'probeId');
clear B;
function [ind, uniGene] = HighExpressionProbes(Genes, Probes, Data);
meanData = mean(Data, 2);
[sGenes, sInd] = sort(Genes);
sProbes = Probes(sInd);
sMean = meanData(sInd);
[uniGene, uniInd] = unique(sGenes);
uniInd = [0; uniInd];
tmpInd = zeros(1, length(uniGene));
for i = 1 : length(uniInd)-1
    [maxV, maxInd] = max(sMean(uniInd(i)+1:uniInd(i+1)));
    tmpInd(i) = uniInd(i) + maxInd;
end;
% [maxV, maxInd] = max(sMean(uniInd(end):end));
% tmpInd(end) = uniInd(end) + maxInd - 1;
ind = sInd(tmpInd);


CMATRIX

%%PAPER

Step 0. l ← 1 where l is the indicator of the levels in the hierarchical system.
w0 ←γmax{w(e): ∀e∈E(G)}whereγ(0<γ<1)is a user specified parameter.

Step 1. (The initial step)
Sort the edge set {e ∈ E(G) : w(e) ≥ w0} as a sequence S = e1,··· ,em such
that w(e1) ≥ w(e2) ≥ ··· ≥ w(em). μ←1,p←0,andLl ←∅.

%%ME

Step 0 + 1. Find the local maximal edges by sorting the set of transformed PCC edge 
values from most similar (maxEdges, largest PCC=shortest edge) to least similar, 
implementing a cutoff threshold to restrict responses to only those most relevant 
to our question.

keepN = 10000;

keepSMALLData = data(sortIndex(1:keepN), :);
keepGeneSym= uniGene(sortIndex(1:keepN));
cMatrix = massivePCC_withoutNaN(keepSMALLData);





MASSIVE PCC (already sorted, so no need after assembled)
% %%%%%%%%%%% Step 1 - Compute PCC matrix %%%%%%%%%%%%%%%%

 % dataMat is a matrix of N-by-M containing NaN entries (N > 1)
% dataVec is a row vector of 1-by-M containing NaN entriesb

% PCC_vec is a column vector of N-by-1 which is the PCC between each row
% of dataMat and dataVec
% minCommonThresh is the minimum number of data pairs required for
% reporting a meaningful PCC value
function PCC_mat = massivePCC_withoutNaN(dataMat);
[nRow, nCol] = size(dataMat);

dataMat(isnan(dataMat)) = 0;


if (nRow > 1)
    
    sumX = sum(dataMat, 2);
    meanX = sumX / nCol;
    sumV = sqrt(sum(dataMat.*dataMat, 2) - (sumX .* sumX)/nCol);
    
    PCC_mat = zeros(nRow, nRow);
    
    for i = 1 : nRow
        if (mod(i, 10000) == 0)
            i
            break
        end;
        PCC_mat(i, i+1:end) = (dataMat(i+1:end,:) * dataMat(i,:)' - sumX(i) * sumX(i+1:end) / nCol)./(sumV(i) * sumV(i+1:end));
    end;
    PCC_mat = PCC_mat + PCC_mat';
    PCC_mat(1 : nRow+1 : end) = 1;
else
    PCC_mat = corr(dataMat(1, ind), dataVec(1, ind));
end;



 

%%%%%%%% Step 2 - identify co-expression modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Algorithm parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = .6; t = 1; lambda = 1;
%%%%%%%% Run the algorithm


LOCAL MAX localMaximumQCM


%%PAPER

Step 2. (Starting a new search). p ← p + 1, Cp ← V (eμ). Ll ← Ll ∪ {Cp}. Step 3. (Grow)
Substep 3.1. If V (G) − V (Cp) = ∅, then go to Step 4, otherwise continue: Pick v ∈ V (G) − V (Cp) such that c(v, Cp) is a maximum.
If
￼c(v, Cp) ≥ αnd(Cp) (1) wheren=|V(Cp)|andαn =1− 1 withλ≥1andt≥1asuserspecified
￼2λ(n+t)
parameters, then Cp ← Cp ∪ {v} and go back to Substep 3.1.
Substep 3.2. μ ← μ + 1. If μ > m go to Step 4. 􏰋p−1
Substep 3.3. Suppose eμ = xy. If at least one of x, y ∈/ i=1 V (Ci) then go to Step 2, otherwise go to Substep 3.2.

%%ME

Step 1. Find the local maximal edges by sorting the set of 
transformed PCC edge values from most similar (maxEdges, largest PCC=shortest edge) 
to least similar, implementing a cutoff threshold to restrict responses to only those 
most relevant to our question.

Step 2. While the clique and subgraph are the same (ie V(G)=V(Cp)), and the length of 
the maximum index (sortMaxInd) is the current search space,

Substep 1    We pick a sorted, maximum vertex (sortMaxV) from the current vertex indices 
(currentIni, t (unique cMatrix coordinate)) of interest that is below the user specified 
gamma threshold multiplied by the max vertex for all overall cliques (gamma* sortMaxV), 
until the graph is more inclusive (by user specified parameters lambda and t) than the 
clique, ie it has added new vertices. This iterative/recursive process dictates the number
 of additional vertices ultimately added to the quasi-clique (V(Cp)) compared to the standard
  clique (V(G)).

Substep 2    Move to next vertex indices (cMatrix coordinate). If the number of vertex 
exceeds that of the clique, proceed to substep 3.

Substep 3    If at least one vertex is not part of the original clique being reviewed, 
ie you have exhausted correlatively expressed genes/nodes, start a new clique evaluation in step 1.

function C = localMaximumQCM(cMatrix, gamma, t, lambda);

C = {};
%%Step 1 - find the local maximal edges
[maxV, maxInd] = max(cMatrix);
maxEdges = []; maxW = [];
for i = 1 : length(maxInd)
    if maxV(i)== max(cMatrix(maxInd(i), :));
        maxEdges = [maxEdges; maxInd(i), i];
        maxW = [maxW; maxV(i)];
    end;
end;

[sortMaxV, sortMaxInd] = sort(maxW, 'descend');
sortMaxEdges = maxEdges(sortMaxInd, :);

length(sortMaxInd)

currentInit = 1; noNewInit = 0;

nodesInCluster = [];

while (currentInit <= length(sortMaxInd)) & (noNewInit == 0)
    if (sortMaxV(currentInit) < gamma * sortMaxV(1))
        noNewInit = 1;
    else
        if (ismember(sortMaxEdges(currentInit, 1), nodesInCluster)==0 & ...
                ismember(sortMaxEdges(currentInit, 2), nodesInCluster)==0)
            newCluster = sortMaxEdges(currentInit, :);
            addingMode = 1;
            currentDensity = sortMaxV(currentInit);
            nCp = 2;
            totalInd = 1 : size(cMatrix, 1);
            remainInd = setdiff(totalInd, newCluster);
            while addingMode == 1
                
                neighborWeights = sum(cMatrix(newCluster, remainInd));
                
                [maxNeighborWeight, maxNeighborInd] = max(neighborWeights);
                c_v = maxNeighborWeight/nCp;
                alphaN = 1 - 1/(2*lambda*(nCp+t));
                if (c_v >= alphaN * currentDensity)
                    newCluster = [newCluster, remainInd(maxNeighborInd)];
                    nCp = nCp+1;
                    currentDensity = (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2);
                    remainInd = setdiff(remainInd, remainInd(maxNeighborInd));
                else
                    addingMode = 0;
                end;
            end;
            nodesInCluster = [nodesInCluster, newCluster];
            C = [C, newCluster];
            
        end;
    end;
    currentInit = currentInit + 1
end;


%%%%%%%% Step 3 - Merge the overlapped networks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%PAPER
Step 4. (Merge)
Substep 4.1.
List all members of Ll as a sequence C1,··· ,Cs such that |V (C1)| ≥ |V (C2)| ≥ · · · ≥ |V (Cs)| where s ← |Ll|.
h ← 2, j ← 1.
Substep 4.2. If |Cj ∩Ch| > βmin(|Cj|,|Ch|) (where β (0 < β < 1) is a user specified parameter), then Cs+1 ← Cj ∪ Ch and the sequence Ll is rearranged as follows
C1,··· ,Cs−1 ← deleting Cj,Ch from C1,··· ,Cs+1 and s ← s − 1, h ← max{h − 2, 1}, and go to Substep 4.4.
Substep 4.3. j ← j + 1. If j < h go to Substep 4.2.
Substep 4.4. h ← h + 1 and j ← 1. If h ≤ s go to Substep 4.2.



%%%ME
Step 3. List all vertexes as a vector sequence, once again in order of decreasing density. 
Sort and iteratively merge networks with overlapping/coexpressed nodes, allowed for by an 
overlap ratio threshold specified by the user by going through all ordered clusters consisting 
of greater than 3 vertex.

Substep 1    If the length between any two vertex is greater than the user specified beta 
threshold, remove the more poorly connected of the two nodes; Move to Step 4.

Substep 2    Evaluate neighboring node’s connections to its surrounding nodes; proceed to 
step 3 substep 1, evaluate length between all node/vertex pairs in.

%%%%%%%% Allowed overlap ratio threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.4; 

%%%%%%%% Sort and iteratively merge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeC = zeros(1, length(C));
for i = 1 : length(C)
    sizeC(i) = length(C{i});
end;

[sortC, sortInd] = sort(sizeC, 'descend');
C = C(sortInd);

ind = find(sortC >= 3);

mergedCluster = C(ind);
mergeOccur = 1; 
currentInd = 0;

while mergeOccur == 1
    mergeOccur = 0;
    while currentInd < length(mergedCluster)
        currentInd = currentInd + 1;
        excludeInd = [];
        if (currentInd < length(mergedCluster))
            keepInd = 1 : currentInd;
            for j = currentInd+1 : length(mergedCluster)
                interCluster = intersect(mergedCluster{currentInd}, mergedCluster{j});
                if length(interCluster) >= beta*min(length(mergedCluster{j}), length(mergedCluster{currentInd}))
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

fid = fopen(strcat(filepath, 'clusterGeneList.txt'), 'w');
for i = 1 : length(mergedCluster)
    fprintf(fid, '%s\t', strcat('cluster ', num2str(i)));
    geneList = keepGeneSym(mergedCluster{i});
    for j = 1 : length(mergedCluster{i})-1
        fprintf(fid, '%s\t', geneList{j+1});
    end;
    fprintf(fid, '\n');
end;
fclose(fid);

%%ME
META GENE T_TEST
%stratified LC patient samples by column number
noResponse= keepSMALLData(:,[13,14,25,16,74,76,40,44,64,68]);

response= keepSMALLData(:,[22,37,17,6,72,94,73,82,8,18]);

%avg expression of strata collapsed by sample
for i=1:1:length(sizeMergedCluster)
temp{i}= response(mergedCluster{i}, :);
%temp{i}=(temp{i})
tempMeanClusterData{i} = mean(temp{i});

temp2{i}= noResponse(mergedCluster{i}, :);
%temp2{i}=(temp2{i})
tempNoMeanClusterData{i} = mean(temp2{i});

%ttest among strata avg/PCA
[h,p] = ttest2 ((tempMeanClusterData{i}),(tempNoMeanClusterData{i}));

%intraclusterPvalue=(p, strcat('cluster' {i} 'p-value', num2str(i))p-value{i+1}));
pv(i)=[p]
end
%keep clusters pv<.05
for i=1:1:length(sizeMergedCluster)
if pv(i)<.1
      intraClusterDiff = [keepGeneSym(mergedCluster{i})]
end


end

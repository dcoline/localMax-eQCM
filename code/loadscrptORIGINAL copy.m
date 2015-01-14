clear;
close all;

filepath = '/Users/danielcolin/Dropbox/OSU/SECOND YEAR/Spring14/Thesis KUN/TWO/code/';
datapath = '/Users/danielcolin/Dropbox/OSU/SECOND YEAR/Spring14/Thesis KUN/TWO/code/';
datafile = 'ccleLCexpr.txt';
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



% load(strcat(filepath, 'matrix.mat'));
% load(strcat(filepath, 'matrixProbeGeneSym.mat'));
% 
% meanData = mean(matrix, 2);
% [sortMean, sortIndex] = sort(meanData, 'descend');
% matrix = matrix(sortIndex, :);
% geneSym = gene(1 : 20067); gene = gene(sortIndex);
% probeId = probe(1: 20067); probe = probe(sortIndex);

[ind, uniGene] = HighExpressionProbes(geneSym, probeId, dataMatrix);
data = dataMatrix(ind, :);
meanData = mean(data, 2);
[sortMean, sortIndex] = sort(meanData, 'descend');


keepN = 10000;

keepSMALLData = data(sortIndex(1:keepN), :);
keepGeneSym= uniGene(sortIndex(1:keepN));
cMatrix = massivePCC_withoutNaN(keepSMALLData);
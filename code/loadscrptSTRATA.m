 clear;
close all;

filepath = '/Users/mor152/Dropbox/OSU/Spring14/Thesis KUN/TWO/code/';
datapath = '/Users/mor152/Dropbox/OSU/Spring14/Thesis KUN/TWO/code/';
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

%stratified LC patient samples by column number

%resp  32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19
response= data(:,[10,59,21,63,81,56,24,26,100,34,55,8,54,82,16,40,16,77,73,17,37,22]);
% non  44,16,3,85,94,6,102,100,82,18,64,68,
noResponse= data(:,[53,112,77,14,73,75,25,108,124,103,149,62,13,94,23,89,64,18,35,47,53,28,44,41,99]);


keepN = 10000;

keepResponse = response(sortIndex(1:keepN), :);
keepNoResponse = noResponse(sortIndex(1:keepN), :);
keepGeneSym= uniGene(sortIndex(1:keepN));
keepSMALLSym = geneSym(1:keepN);
keepSMALLProbe = probeId(1:keepN);
 
cMatrix = massivePCC_withoutNaN(keepResponse);
cMatrix2 = massivePCC_withoutNaN(keepNoResponse);
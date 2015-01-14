 clear;
close all;

filepath = '/Users/Daniel/Downloads/TWO/code/';
datapath = '/Users/Daniel/Downloads/TWO/code/';
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
%response= data([11816,11814,11815,11812,11811,1810],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
%response= data([7351,7120,7345,7346,7348,7349],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
response= data([7741,4382,7739,7740,4383],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
%noResponse= data([11820,11818,11819,11816,11815,1814],[13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
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
MT= data([11816,11814,11815,11812,11811,1810],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
HIST= data([7351,7120,7345,7346,7348,7349],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
IFIT= data([7749,7750,4383,7742,7740,7741,4384],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
LOC= data([9440,10514,12526],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
ATA= data([1110,1111],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
NFIA= data([12308,8442],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
%IFITM= data([7749,7750],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
RANDOM= data([15961,10740,3867,15093,14042,6757,1581,10848,10451,17216],[32,74,57,42,53,63,12,59,43,15,25,81,95,59,28,19,21,63,55,8,54,82,6,40,16,77,73,17,37,22,13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);

%noResponse= data([11820,11818,11819,11816,11815,1814],[13,72,44,16,3,85,94,6,102,100,82,18,64,68,76,14,72,74,25,102,13]);
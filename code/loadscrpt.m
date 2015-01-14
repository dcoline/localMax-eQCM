clear;
close all;

filepath = '/Users/danielcolin/Dropbox/OSU/SECOND YEAR/Spring14/Thesis KUN/TWO/code/';
datapath = '/Users/danielcolin/Dropbox/OSU/SECOND YEAR/Spring14/Thesis KUN/TWO/code/';
%datafile = 'ccleLCexpr.txt';
datafile = 'CCLE_EXP_FILT.txt';
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

cMatrix = massivePCC_withoutNaN(dataMatrix);
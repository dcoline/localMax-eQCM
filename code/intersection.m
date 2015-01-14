clear;
close all;
%directory
filepath = '/Users/mor152/Dropbox/OSU/Spring14/Thesis KUN/TWO/code/';
datapath = '/Users/mor152/Dropbox/OSU/Spring14/Thesis KUN/TWO/code/';

%read nonresponsive clusters
[num,raw]=xlsread('clusterGeneList11.xlsx');
noResp=raw;

%each column is cluster
for i=1:1:length(noResp);
     clusterNo{i}=noResp(:,i);
end;

%read responsive clusters
[num,raw]=xlsread('clusterGeneList22.xlsx');
resp=raw;

%each column is cluster
for i=1:1:length(resp);
     clusterResp{i}=resp(:,i); %first strata must be larger in both X&Y
end;




%for clusters >3, check overlap btwn clusters
for j=(1:123);
    for k=(1:298);
    B{j,k}=intersect(clusterResp{j},clusterNo{k});    
    C{j,k}=size(intersect(clusterResp{j},clusterNo{k}));
    D{j,k}=union(clusterResp{j},clusterNo{k});
    E{j,k}=size(union(clusterResp{j},clusterNo{k}));
   
    end;
end;

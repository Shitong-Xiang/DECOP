clear;clc

MID_Perf = readtable('/share/inspurStorage/home1/ISTBI_data/ABCD_MID/Behaviour_Score/ABCD_MID_Perf_Site.txt');
load('/share/inspurStorage/home1/ISTBI_data/ABCD_MID/Behaviour_Score/Covas.mat');
FD_EX = find(Cova_Site.MeanFD > 0.5);
Cova_Site(FD_EX,:) = [];
MEs = dir('/share/inspurStorage/home1/ISTBI_data/ABCD_MID/Exact_vs_Relative_Value/*MEs.txt');
name = MID_Perf.Properties.VariableNames;

Salience_ind = readtable('/share/inspurStorage/home1/ISTBI_data/ABCD_MID/Exact_vs_Relative_Value/Visual_Value/Salience.csv');
Salience_ind = find(Salience_ind.Ratio > 0);
Value_ind = readtable('/share/inspurStorage/home1/ISTBI_data/ABCD_MID/Exact_vs_Relative_Value/Visual_Value/Value.csv');
Value_ind = find(Value_ind.Ratio > 0);

[~,ind1,ind2] = intersect(Cova_Site.SubID,MID_Perf.SubID);
Cova_reg = [ones(length(ind1),1),table2array(Cova_Site(ind1,2:end))];
[~,~,Y] = regress(MID_Perf.MeanRT(ind2),Cova_reg);

tmp_MEs = readtable(fullfile(MEs(1).folder,MEs(1).name));
% tmp_MEs = readtable(fullfile(MEs(2).folder,MEs(2).name));
MEs_names = tmp_MEs.Properties.VariableNames;
MEs_names1 =  MEs_names(Salience_ind);
tmp_MEs1 = tmp_MEs(:,Salience_ind);
% MEs_names1 =  MEs_names(Value_ind);
% tmp_MEs1 = tmp_MEs(:,Value_ind);
tmp_MEs1 = mat2cell(table2array(tmp_MEs1(ind1,:)),length(ind1)*ones(1,1),ones(1,size(tmp_MEs1,2)));
[~,~,X] = cellfun(@(x) regress(x,Cova_reg),tmp_MEs1,'un',0);
X1 = cell2mat(X);
tmp_MEs2 = tmp_MEs;
MEs_names2 = MEs_names;
tmp_MEs2(:,Salience_ind) = [];
MEs_names2(Salience_ind) = [];
% tmp_MEs2(:,Value_ind) = [];
% MEs_names2(Value_ind) = [];
tmp_MEs2 = mat2cell(table2array(tmp_MEs2(ind1,:)),length(ind1)*ones(1,1),ones(1,size(tmp_MEs2,2)));
[~,~,X] = cellfun(@(x) regress(x,Cova_reg),tmp_MEs2,'un',0);
X2 = cell2mat(X);

names1 = ['Perf',MEs_names1];
names2 = ['Perf',MEs_names2];
obs_diff = x_R_square_diff(Y,X1,X2,names1,names2);

[~,bootsam]  = bootstrp(10000,@(x) x,1:length(ind1));
bootstat = zeros(1,10000);
tic
parfor i = 1:10000
    disp(i);
    bootstat(i) = x_R_square_diff(Y(bootsam(:,i)),X1(bootsam(:,i),:),X2(bootsam(:,i),:),names1,names2);
end
toc

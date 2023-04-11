clear;clc

load('H:\ABCD\ABCD_MID\CCA_analysis\Covas.mat');
Cova_Site(Cova_Site.MeanFD > 0.5,:) = [];

Value = readtable('H:\ABCD\ABCD_MID\CCA_analysis\Value_pos_WVCNA_MEs.txt');
Salience = readtable('H:\ABCD\ABCD_MID\CCA_analysis\Salience_pos_WVCNA_MEs.txt');

Pheno = readtable('H:\ABCD\ABCD_MID\CCA_analysis\All_Domains.csv');
[~,ind1,ind2] = intersect(Cova_Site.SubID,Pheno.SubID);

Cova_Site = Cova_Site(ind1,:);
Value = Value(ind1,:);
Salience = Salience(ind1,:);
Pheno = Pheno(ind2,:);

Cova_reg = [ones(size(Cova_Site,1),1),table2array(Cova_Site(:,2:6))];
Value_MEs = mat2cell(table2array(Value),size(Cova_Site,1)*ones(1,1),ones(1,size(Value,2)));
[~,~,Value_Reg] = cellfun(@(x) regress(x,Cova_reg),Value_MEs,'un',0);
Value_Reg = cell2mat(Value_Reg);
Value_Reg(isnan(Value_Reg)) = 0;
Salience_MEs = mat2cell(table2array(Salience),size(Cova_Site,1)*ones(1,1),ones(1,size(Salience,2)));
[~,~,Salience_Reg] = cellfun(@(x) regress(x,Cova_reg),Salience_MEs,'un',0);
Salience_Reg = cell2mat(Salience_Reg);
Salience_Reg(isnan(Salience_Reg)) = 0;

Beh = mat2cell(table2array(Pheno(:,2:end)),size(Cova_Site,1)*ones(1,1),ones(1,size(Pheno,2)-1));
[~,~,Beh_Reg] = cellfun(@(x) regress(x,Cova_reg),Beh,'un',0);
Beh_Reg = cell2mat(Beh_Reg);

[A1,B1,r1,U1,V1,stats1] = canoncorr(Value_Reg,Beh_Reg(:,[11:12,22:23,26:end]));
[A2,B2,r2,U2,V2,stats2] = canoncorr(Salience_Reg,Beh_Reg(:,[11:12,22:23,26:end]));
Eta1 = 1- exp(sum(log(1-power(r1,2))));
Eta2 = 1- exp(sum(log(1-power(r2,2))));

parfor i = 1:10000
    disp(i);
    perm_ind = randperm(size(Beh_Reg,1));
    [~,~,r01(i,:)] = canoncorr(Value_Reg,Beh_Reg(perm_ind,[11:12,22:23,26:end]));
    Eta01(i) = 1- exp(sum(log(1-power(r01(i,:),2))));
    [~,~,r02(i,:)] = canoncorr(Salience_Reg,Beh_Reg(perm_ind,[11:12,22:23,26:end]));
    Eta02(i) = 1- exp(sum(log(1-power(r02(i,:),2))));
end

clear;clc

load('H:\ABCD\ABCD_MID\CCA_analysis\CCA_Output.mat');

Stats(1,1) = 1-(1-Eta1)/(1-mean(Eta01));
Stats(1,2) = sum(Eta01 > Eta1)/10000;
Stats(1,3) = 1-(1-Eta2)/(1-mean(Eta02));
Stats(1,4) = sum(Eta02 > Eta2)/10000;

Stats(2:9,1) = 1-(1-power(r1,2))./(1-power(mean(r01),2));
Stats(2:9,2) = sum(r01 > r1)/10000;
Stats(2:9,3) = 1-(1-power(r2,2))./(1-power(mean(r02),2));
Stats(2:9,4) = sum(r02 > r2)/10000;
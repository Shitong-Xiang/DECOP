clear;clc

Value = '/share/inspurStorage/home1/ISTBI_data/ABCD_V2_MID/Experimental_decode';
Salience = '/share/inspurStorage/home1/ISTBI_data/ABCD_V2_MID/Experimental_decode';
sub = dir(fullfile(Value,'sub-*'));

load('/share/inspurStorage/home1/ISTBI_data/ABCD_V2_MID/Experimental_decode/x_Simulation_TS.mat');

tic
for j = 1:1000
    select = randi(length(TS_R),1000,1);
    Value_sim_beta = 0.5+randn(1000,1);
    Salience_sim_beta = 0.5+randn(1000,1);
    Noise = randn(1000,822);
    disp(j);
    parfor i = 1:1000
        disp(i);
        sub_TS = Value_sim_beta(i)*Value_TS(select(i),:)+Salience_sim_beta(i)*Salience_TS(select(i),:)+Noise(i,:);
        b = glmfit(Value_TS(select(i),:)',sub_TS');
        uni1(i,1) = b(2);
        b = glmfit(Salience_TS(select(i),:)',sub_TS');
        uni2(i,1) = b(2);
    end
    r(j) = corr(uni1,uni2);
end
toc
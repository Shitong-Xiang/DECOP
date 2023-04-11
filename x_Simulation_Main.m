clear;clc
base = [0 0 0 1 2];
step = -1:0.1:1;
x1 = [2 1 0 0 0];
x2 = [0 0 0 1 2];
X = [x1',x2'];
parfor time = 1:1000
    noise = randn(1000,5);
    [r_1(:,time),L_t_1(:,time),V_t_1(:,time),r_2(:,time),L_t_2(:,time),V_t_2(:,time),r_3(:,time),L_t_3(:,time),V_t_3(:,time),r_4(:,time),L_t_4(:,time),V_t_4(:,time)]= x_one_simulation(step,noise,X);
end

stat(:,1) = mean(r_1,2);
stat(:,2) = mean(L_t_1,2);
stat(:,3) = mean(V_t_1,2);
stat(:,4) = mean(r_2,2);
stat(:,5) = mean(L_t_2,2);
stat(:,6) = mean(V_t_2,2);
stat(:,7) = mean(r_3,2);
stat(:,8) = mean(L_t_3,2);
stat(:,9) = mean(V_t_3,2);
stat(:,10) = mean(r_4,2);
stat(:,11) = mean(L_t_4,2);
stat(:,12) = mean(V_t_4,2);

sig = sig';
function [r_1,L_t_1,V_t_1,r_2,L_t_2,V_t_2,r_3,L_t_3,V_t_3,r_4,L_t_4,V_t_4]= x_one_simulation(step,noise,X)

for j = 1:length(step)
    disp(j);
    sig{j} = round([2*step(j),step(j),0,0,0]+[0 0 0 1 2],2);
    a1 = 1+randn(1000,1);
    y1 = a1*sig{j}+noise;
    a2 = 0.05*(j-1)+randn(1000,1);
    a3 = 1-0.05*(j-1)+randn(1000,1);
    y2 = a2*[2 1 0 1 2]+a3*[-2 -1 0 1 2]+noise;
    a4 = 0.05*(j-1) + randn(1000,1);
    a5 = a4+randn(1000,1);
    a6 = 1-a4+randn(1000,1);
    y3 = a5*[2 1 0 1 2]+a6*[-2 -1 0 1 2]+noise;
    a7 = -1+0.1*(j-1)+randn(1000,1);
    a8 = 1+randn(1000,1);
    y4 = a7*[2 1 0 0 0]+a8*[0 0 0 1 2]+noise;
    for i = 1:1000
        b1(:,i) = glmfit(X,y1(i,:)');
    end
    r_1(j,1) = corr(b1(2,:)',b1(3,:)');
    [h p c t] = ttest(b1(2,:)');
    L_t_1(j,1) = t.tstat;
    [h p c t] = ttest(b1(3,:)');
    V_t_1(j,1) = t.tstat;
    for i = 1:1000
        b3(:,i) = glmfit(X,y3(i,:)');
    end
    r_2(j,1) = corr(b3(2,:)',b3(3,:)');
    [h p c t] = ttest(b3(2,:)');
    L_t_2(j,1) = t.tstat;
    [h p c t] = ttest(b3(3,:)');
    V_t_2(j,1) = t.tstat;
    for i = 1:1000
        b2(:,i) = glmfit(X,y2(i,:)');
    end
    r_3(j,1) = corr(b2(2,:)',b2(3,:)');
    [h p c t] = ttest(b2(2,:)');
    L_t_3(j,1) = t.tstat;
    [h p c t] = ttest(b2(3,:)');
    V_t_3(j,1) = t.tstat;
    for i = 1:1000
        b4(:,i) = glmfit(X,y4(i,:)');
    end
    r_4(j,1) = corr(b4(2,:)',b4(3,:)');
    [h p c t] = ttest(b4(2,:)');
    L_t_4(j,1) = t.tstat;
    [h p c t] = ttest(b4(3,:)');
    V_t_4(j,1) = t.tstat;
end
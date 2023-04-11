function diff =  x_R_square_diff(Y,X1,X2,name1,name2)

glm_table = array2table([Y,X1],'VariableNames',name1);
model_set = [name1{1},'~',name1{2}];
for k = 3:length(name1)
    model_set = [model_set,'+',name1{k}];
end
glm = fitlm(glm_table,model_set);
R1 = glm.Rsquared.Adjusted;

glm_table = array2table([Y,X2],'VariableNames',name2);
model_set = [name2{1},'~',name2{2}];
for k = 3:length(name2)
    model_set = [model_set,'+',name2{k}];
end
glm = fitlm(glm_table,model_set);
R2 = glm.Rsquared.Adjusted;
diff = R1 - R2;


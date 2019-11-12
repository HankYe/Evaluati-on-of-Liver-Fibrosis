function [T1,T2,dec_value,TP, FP, TN, FN]=myclassifysvmROC_CV_libsvm(X_test,T_teacher,X_train,T_train,C_value, g_value)


% SVM train
cmd = [ '-t 2', ' -c ',num2str(C_value), ' -g ', num2str(g_value), ' -b 1'];
model = svmtrain(T_train, X_train,cmd); %PCBP 80

% SVM predict value
[output, acc, dec] = svmpredict(T_teacher, X_test,model,'-b 1');
output_test = output;
dec_value = dec(:,1);


% 计算分类结果
T1 = output_test; % test result
T1(find(output_test==-1))=0; % benign re-label 0
T2 = T_teacher; % real result
T2(find(T_teacher==-1))=0; % benign re-label 0

% 
TP = sum(T1&T2); % real - m, test - m
FP = sum(T1&~T2);% real - b, test - m
FN = sum(~T1&T2);% real - m, test - b
TN = sum(~T1&~T2); % real - b, test - b  





end


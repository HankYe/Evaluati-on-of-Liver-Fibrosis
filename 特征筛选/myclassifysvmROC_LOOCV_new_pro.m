%function [sens_predict, spes_predict, accuracy,sensitivity,specificity,PPV,NPV,MCC ] = myclassifysvmROC_LOOCV_new_pro( X_train,X_test,T,C_value,g_value)
function [sens_predict, spes_predict,mean_accu,mean_sens,mean_spes,mean_ppv,mean_npv,mean_mcc]= myclassifysvmROC_LOOCV_new_pro( X_train,X_test,T,C_value,g_value)
% Leave-one-out cross validation(LOOCV) of svm classification (留一法)
% X_train: 训练集
% X_test: 测试集
% T： 分类目标
% C_value: C for SVM
% g_value: gamma for SVM

e=1e-10;
sampletotal=size(X_train,1);

M = find(T==1); % malignant index 
B = find(T==0); % benign index 


n_malignant = length(M); % malignant
n_benign = sampletotal- n_malignant; % benign

% Rearrange X, [69 malignant, 69 benign]
% Mi= M(randperm(n_malignant));
X_train_scale(1:n_malignant,:) = X_train(M,:);
X_test_scale(1:n_malignant,:) = X_test(M,:);

% Bi= B(randperm(n_benign));
X_train_scale(n_malignant+1:sampletotal,:) = X_train(B,:);
X_test_scale(n_malignant+1:sampletotal,:) = X_test(B,:);


% Rearranged Target vector, 1 for malignant, -1 for benign
T_scale=[ones(1,n_malignant) -ones(1,n_benign)]'; % [malignant benign] 138X1
output_test = zeros(sampletotal,1); % 138X1
dec_value = zeros(sampletotal,1); % 138X1
% acc = zeros(sampletotal,1); % 138X1

for k = 1:sampletotal
    
    if k == 1
        D_train = X_train_scale(2:sampletotal,:);
        T_train = T_scale(2:sampletotal);
    elseif k == sampletotal
        D_train = X_train_scale(1:sampletotal-1,:);
        T_train = T_scale(1:sampletotal-1);
    else
        D_train = [X_train_scale(1:k-1,:)' X_train_scale(k+1:sampletotal,:)']';
        T_train = [T_scale(1:k-1)' T_scale(k+1:sampletotal)']';
    end
    
    D_test = X_test_scale(k,:);
    T_test = T_scale(k);
    % SVM train
    cmd = [ '-t 2', ' -c ',num2str(C_value), ' -g ', num2str(g_value), ' -b 1'];
    model = svmtrain(T_train, D_train,cmd); %PCBP 80
%     model = svmtrain(T_train, D_train, '-t 2 -c 5.6569 -g 0.0039063 -b 1'); %PCBP 144
%     model = svmtrain(T_train, D_train, '-t 2 -c 2 -g 0.0884 -b 1'); %PCBP 80
%     model = svmtrain(T_train, D_train, '-t 2 -c 2.8284 -g 0.0313 -b 1'); %PCBP 80
%     model = svmtrain(T_train, D_train, '-t 2 -c 724.0773 -g 0.0110485 -b 1'); %Ranklet144
%     model = svmtrain(T_train, D_train, '-t 2 -c 1024 -g 0.0110 -b 1'); %Ranklet72
%     model = svmtrain(T_train, D_train, '-t 2 -c 0.5000 -g 0.0442 -b 1'); % GLCM
%     model = svmtrain(T_train, D_train, '-t 2 -c 2.8284 -g 1 -b 1'); %LBP
%     model = svmtrain(scale_label_train, scale_inst_train, '-c 1 -g 1 -b 1');

%     model = svmtrain(T_train, D_train, '-t 2 -c 0.7071 -g 0.3536 -b 1'); %PCBP o=4 s=3
%     model = svmtrain(T_train, D_train, '-t 2 -c 0.5000 -g 0.3536 -b 1'); %PCBP o=4 s=4
%     model = svmtrain(T_train, D_train, '-t 2 -c 5.6569 -g 0.0625 -b 1'); %PCBP o=4 s=5
%     model = svmtrain(T_train, D_train, '-t 2 -c 2.8284 -g 0.0884 -b 1'); %PCBP o=4 s=6

%     model = svmtrain(T_train, D_train, '-t 2 -c 0.7071 -g 0.0884 -b 1'); %PCBP o=6 s=3
%     model = svmtrain(T_train, D_train, '-t 2 -c 32 -g 0.0313 -b 1'); %PCBP o=6 s=4
%     model = svmtrain(T_train, D_train, '-t 2 -c 22.6274 -g 0.0156 -b 1'); %PCBP o=6 s=5
%     model = svmtrain(T_train, D_train, '-t 2 -c 11.3137 -g 0.0221 -b 1'); %PCBP o=6 s=6

%     model = svmtrain(T_train, D_train, '-t 2 -c 8 -g 0.0884 -b 1'); %PCBP o=8 s=3
%     model = svmtrain(T_train, D_train, '-t 2 -c 1 -g 0.1768 -b 1'); %PCBP o=8 s=4
%     model = svmtrain(T_train, D_train, '-t 2 -c 0.7071 -g 0.0221 -b 1'); %PCBP o=8 s=5

%     model = svmtrain(T_train, D_train, '-t 2 -c 2 -g 0.0884 -b 1'); %PCBP o=10 s=3
%     model = svmtrain(T_train, D_train, '-t 2 -c 2 -g 0.1250 -b 1'); %PCBP o=10 s=4
%     model = svmtrain(T_train, D_train, '-t 2 -c 1 -g 0.0625 -b 1'); %PCBP o=10 s=5
%     model = svmtrain(T_train, D_train, '-t 2 -c 1.4142 -g 0.0884 -b 1'); %PCBP o=10 s=6
    
%     model = svmtrain(T_train, D_train, '-t 2 -c 1 -g 0.1250 -b 1'); %PCBP o=12 s=3
%     model = svmtrain(T_train, D_train, '-t 2 -c 1.4142 -g 0.0884 -b 1'); %PCBP o=12 s=4
%     model = svmtrain(T_train, D_train, '-t 2 -c 0.25 -g 0.0884 -b 1'); %PCBP o=12 s=5
%     model = svmtrain(T_train, D_train, '-t 2 -c 90.5097 -g 0.0028 -b 1'); %PCBP o=12 s=6
    [output, acc, dec] = svmpredict(T_test, D_test,model,'-b 1');
    output_test(k) = output;
    dec_value(k) = dec(:,1);
end

%  ROC calculation
  dec_value_sort = sort(dec_value);
for time = 1:length(dec_value)
    thresh0 = dec_value_sort(time);
    for time00 = 1:length(dec_value)
        if (dec_value(time00)-thresh0>0)
            predict_label(time,time00)=1;
        else
            predict_label(time,time00)=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = predict_label;
% t2 = T_teacher;
% t2(find(T_teacher==-1))=0;
t2 = [ones(1,n_malignant) zeros(1,n_benign)]; % real result;

sens_predict = zeros(1,length(dec_value));
spes_predict = zeros(1,length(dec_value));
for time = 1:length(dec_value)
    TP(time) = sum(t1(time,:)&t2);
    FP(time) = sum(t1(time,:)&~t2);
    FN(time) = sum(~t1(time,:)&t2);
    TN(time) = sum(~t1(time,:)&~t2);
    sens_predict(time) = TP(time)/(TP(time)+FN(time));
    spes_predict(time) = TN(time)/(TN(time)+FP(time));
    accu(time)=(TP(time)+TN(time))/(TP(time)+FP(time)+FN(time)+TN(time));
    sens(time)=TP(time)/(TP(time)+FN(time));
    spes(time)=TN(time)/(TN(time)+FP(time));
    ppv(time)=TP(time)/(TP(time)+FP(time));
    npv(time)=TN(time)/(TN(time)+FN(time));
    mcc(time)=(TP(time)*TN(time)-FP(time)*FN(time))/sqrt((TP(time)+...
        FP(time)) *(TP(time)+FN(time))*(TN(time)+FP(time))*(TN(time)+FN(time)) );
    
end

mean_accu=mean(accu);
mean_sens=mean(sens);
mean_spes=mean(spes);
mean_ppv=mean(ppv);
mean_npv=mean(npv);
mean_mcc=mean(mcc);

std_accu=std(accu);
std_sens=std(sens);
std_spes=std(spes);
std_ppv=std(ppv);
std_npv=std(npv);
std_mcc=std(mcc);


%原蔡凌云程序
%     T1 = output_test; % test result
%     T1(find(output_test==-1))=0; % benign re-label 0
%     T2 = [ones(1,n_malignant) zeros(1,n_benign)]'; % real result
% 
%  
%     TP = sum(T1&T2); % real - m, test - m
%     FP = sum(T1&~T2);% real - b, test - m
%     FN = sum(~T1&T2);% real - m, test - b
%     TN = sum(~T1&~T2); % real - b, test - b  
% 
%     accuracy = (TP+TN)/(TP+TN+FP+FN+e);
%     sensitivity = TP/(TP+FN+e);
%     specificity = TN/(TN+FP+e);
%     PPV = TP/(TP+FP+e);
%     NPV = TN/(TN+FN+e);
%     MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)+e);

end


function [mean_ACC] = libsvm_CV(feature_input, label_target, K_fold,temp,g_sel,C_sel)
% load Target138.mat;
% load Breast138_morph11_auto_new.mat;
% load PC_LBPV138_8ori_6scale_r2p16.mat
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load T34-1.mat
% load Feature34.mat

% load GAresult_gen30_fun1_4; %GA_mRMR结果
% load GAresult_gen30_fun0_1; %GA结果

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 x = [feature_input'];
% x = [Breast138_morph11_auto_new' ; PC_LBPV138_8ori_6scale_r2p16'];
% x = [PC_LBPV138_8ori_6scale_r2p16'];
% x = [Breast138_morph11_auto_new'];

[no dim] = size(x');  % database

% feature normalization [-1,1]
for i=1:dim
    x(i,:)=((x(i,:)-min(x(i,:)))/(max(x(i,:))-min(x(i,:))))*2-1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = x';
T = label_target;
% T = T138;
K_fold = K_fold;

sampletotal=size(X,1);

M = find(T==1); % malignant index 
B = find(T==0); % benign index 


n_malignant = length(M); % malignant
n_benign = sampletotal- n_malignant; % benign

% Rearrange X, [69 malignant, 69 benign]
Mi= M(randperm(n_malignant));
X_malignant = X(Mi,:);
T_malignant = ones(1,n_malignant)';
%%% X_scale(1:n_malignant,:) = X_malignant;
%%% %%%此句本来是有的，但是好像没有用到，删了之后调试通过，所以就删了
 
Bi= B(randperm(n_benign));
X_benign = X(Bi,:);
T_benign = -ones(1,n_benign)';

fea_sel = find(temp==1); %%%特征选择所选择特征的序号
% fea_sel = [1:dim];
X_malignant_scale = X_malignant(:,fea_sel);
X_benign_scale = X_benign(:,fea_sel);

% % default setting
% C_sel = 1;
% g_sel = 1/dim; % 1/#features

% indices_m = crossvalind('Kfold',n_malignant,K_fold);
% indices_b = crossvalind('Kfold',n_benign,K_fold);
for i = 1:100
indices_m = crossvalind('Kfold',n_malignant,K_fold);
indices_b = crossvalind('Kfold',n_benign,K_fold);

for m = 1:K_fold
%     test = (indices == m);
%     traininput = X(~test,:);
%     traintarget = T(~test);
%     testinput = X(test, :);
%     testtarget = T(test);
    test_m = (indices_m == m);
    traininput_m = X_malignant_scale(~test_m,:);
    traintarget_m = T_malignant(~test_m);
    testinput_m = X_malignant_scale(test_m, :);
    testtarget_m = T_malignant(test_m);
    
    test_b = (indices_b == m);
    traininput_b = X_benign_scale(~test_b,:);
    traintarget_b = T_benign(~test_b);
    testinput_b = X_benign_scale(test_b, :);
    testtarget_b = T_benign(test_b);
    
    testinput = [testinput_m' testinput_b']';
    testtarget = [testtarget_m' testtarget_b']';
    traininput = [traininput_m' traininput_b']';
    traintarget = [traintarget_m' traintarget_b']';
    
    fold_size(m) = size(testtarget,1);
    [T1,T2,dec_V, tp, fp, tn, fn]=myclassifysvmROC_CV_libsvm(testinput,testtarget,traininput,traintarget,C_sel,g_sel);
    TP_temp(m) = tp;
    FP_temp(m) = fp;
    TN_temp(m) = tn;
    FN_temp(m) = fn;
    

%     sensROC = sens_predict;
%     spesROC = spes_predict;
%     AUC(i,m) = -sum(0.5*(sensROC(2:end)+sensROC(1:end-1)).*(spesROC(2:end) - spesROC(1:end-1)));
    if m ==1
        dec_value(1:fold_size(m)) = dec_V;
        target_value(1:fold_size(m)) = testtarget;
    else
        dec_value(sum(fold_size(1:m-1))+1:sum(fold_size(1:m))) = dec_V;
        target_value(sum(fold_size(1:m-1))+1:sum(fold_size(1:m))) = testtarget;
    end
    
end

TP(i) = sum(TP_temp);
FP(i) = sum(FP_temp);
TN(i) = sum(TN_temp);
FN(i) = sum(FN_temp);

ACC(i) = (TP(i)+TN(i))/(TP(i)+TN(i)+FP(i)+FN(i));
SENS(i) = TP(i)/(TP(i)+FN(i));
SPEC(i) = TN(i)/(TN(i)+FP(i));
PPV(i) = TP(i)/(TP(i)+FP(i));
NPV(i) = TN(i)/(TN(i)+FN(i));

MCC(i) = (TP(i)*TN(i)-FP(i)*FN(i))/sqrt((TP(i)+FP(i))*(TP(i)+FN(i))*(TN(i)+FP(i))*(TN(i)+FN(i)));

%计算AUC K-fold CV
dec_value_sort = sort(dec_value);
for time0 = 1:length(dec_value)
    thresh0 = dec_value_sort(time0);
    for time00 = 1:length(dec_value)
        if (dec_value(time00)-thresh0>0)
            predict_label(time0,time00)=1;
        else
            predict_label(time0,time00)=0;
        end
    end
end

%% ROC
t1 = predict_label;
t2 = target_value; %%%顺序换了？？？？
t2(find(target_value==-1))=0;

sens_predict = zeros(1,length(dec_value));
spes_predict = zeros(1,length(dec_value));
for time0 = 1:length(dec_value)
    TP_predict(time0) = sum(t1(time0,:)&t2);
    FP_predict(time0) = sum(t1(time0,:)&~t2);
    FN_predict(time0) = sum(~t1(time0,:)&t2);
    TN_predict(time0) = sum(~t1(time0,:)&~t2);
    sens_predict(time0) = TP_predict(time0)/(TP_predict(time0)+FN_predict(time0));
    spes_predict(time0) = TN_predict(time0)/(TN_predict(time0)+FP_predict(time0));
end

sensROC(i,:) = sens_predict;
spesROC(i,:) = spes_predict;
v1 = sens_predict; v2 = 1-spes_predict;
AUC(i) = -sum(0.5*(v1(2:end)+v1(1:end-1)).*(v2(2:end) - v2(1:end-1)));

end

mean_ACC=mean(ACC);

 
sensSVM=mean(sensROC,1);
spesSVM=mean(spesROC,1);

plot(1-spesSVM,sensSVM,'LineWidth',0.6,'Color','b');
xlabel('False Positive Fraction'); ylabel('True Positive Fraction');

AUC = -sum(0.5*(sensSVM(2:end)+sensSVM(1:end-1)).*(spesSVM(2:end) - spesSVM(1:end-1)));
fprintf('\nAUC: %g',AUC);
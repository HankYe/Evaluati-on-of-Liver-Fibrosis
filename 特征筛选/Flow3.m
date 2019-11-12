clear all
clc

%% load data
% load('TI_Features_324p2GAmrmr_.mat')
% load('TI_Label_169.mat')
% load('Feature_TI_System_521_0628.mat')
% load('Label_TI_System_521.mat')
load('yesornoforsample_see521to299.mat')
load('yesornoforsample_521to430.mat')
load('yesornoforsample_506to181.mat')
load('YesOrNoForPatients.mat')
load('yesornoforsample_506to181.mat')
load('yesornoforall_1142.mat')
load('Feature_TI_System_521_0626.mat', 'Feature_TI_System_521_0626')
load('Feature_TI_System_521_0626.mat', 'Label_TI_System_521')
load('matlabdata20160710mat.mat', 'index_313to637_of181')
load('matlabdata20160710mat.mat', 'index_313to637to1142_of181')

%% set
all_feature=Feature_TI_System_521_0626(:,find(yesornoforfeature==1));
all_label=Label_TI_System_521+1 ; %0,1
method_p='ttest2'; %choose from  'ttest2'  'ranksum'
% method_featureselection='SRC'; % SRC SRCnew LASSO mrmr Fisher_score reliefF F_score
select_feature_number=40;
k_bootstrap=100; %500;
iteration=3;
flag=1;
p_threshold=0.05;
% %% scale 归一化
% [n_x,n_y]=size(feature);
% 
% % Feature normalization [-1,1]
% for i=1:n_y
%     input_scale(:,i)=((feature(:,i)-min(feature(:,i)))/(max(feature(:,i))-min(feature(:,i))))*2-1;
% end

for k=1:k_bootstrap
%% 数据集构造
% 特征筛选训练集/训练集  9/10
% 测试集  1/10
% 验证集 
clear Train_feature_scale
clear Train_label
clear Test_feature
clear Test_label
clear Train_feature_scale
clear Test_feature_scale
clear output_test
clear P
clear H2
clear index
clear featureselection_index

% [A,B,Train_feature,Train_label,Test_feature,Test_label]=ten_flod(all_feature,all_label);
%label 1 2
[all_feature_scale,A0,B0] = scaling(all_feature);
% Test_feature_scale = scaling(Test_feature,1,A0,B0); %归一化

% %每次bootstrap十折得到的数据都存为一层数据，存一个整体的数据
% All_Train_feature_scale(:,:,k)=Train_feature_scale;
% All_Train_label(:,:,k)=Train_label;
% All_Test_feature_scale(:,:,k)=Test_feature_scale;
% All_Test_label(:,:,k)=Test_label;


all_feature2=[all_feature(find(all_label==1),:)' all_feature(find(all_label==2),:)']';
L_N=length(  find(all_label==1) ) ;
L_P=length(  find(all_label==2) );
all_label2=[ ones(1, L_N  )  1+ones(1,L_P)  ]';  

%% 9/10的数据    特征筛选训练集
[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,all_feature2,all_label2); %4改好的svm但是有nan  3 libsvm
Result_origin(1,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]

%% 单因素分析 p值
[output_test,P,H2,index_p]=p_value_function(Train_feature_scale,Train_label,p_threshold,method_p);
%p_index;
output_test=Train_feature_scale(:,index_p);

All_P(k,:)=P;

[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,output_test,Train_label); %4改好的svm但是有nan  3 libsvm
Result_p(1,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]

%% 特征选择1 
[temp_Feature_selection,featureselection_index_p2GA_mRMR,Data_sort]=feature_selection(output_test,Train_label,'GA_mRMR');
% method= GA_mRMR    SRCnew  SRC  LASSO  mrmr  Fisher_score  reliefF  F_score

%feature_selection_index对应初始特征集位置的特征索引
% GA_mRMR_index=(~H2)*size(all_feature,2);
index_GA_mRMR=zeros(1,size(featureselection_index_p2GA_mRMR,2));
for i=1:size(featureselection_index_p2GA_mRMR,2)
    index_GA_mRMR(i)=index_p(featureselection_index_p2GA_mRMR(i));
end

gamrmr_output=Train_feature_scale(:,index_GA_mRMR ) ;
All_GA_mRMR_index(k,:)=index_GA_mRMR;

[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,gamrmr_output,Train_label); %4改好的svm但是有nan  3 libsvm
Result_gamrmr_output(1,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]

%% 特征选择2
[~,filp_index,Data_sort]=feature_selection(gamrmr_output,Train_label,'SRCnew');
% method= GA_mRMR    SRCnew  SRC  LASSO  mrmr  Fisher_score  reliefF  F_score

%feature_selection_index对应初始特征集位置的特征索引
index_feature_selection=zeros(1,size(filp_index,2));
for i=1:size(filp_index,2)
    index_feature_selection(i)=index_GA_mRMR(filp_index(i));
end


Data_sort=Train_feature_scale(:,index_feature_selection ) ;
% Temp_Feature_selection(k,:)= temp_Feature_selection;
% All_index(k,:)=featureselection_index;
clear All_index_feature_selection;
All_index_feature_selection(k,:)=index_feature_selection;

for i=5:5:50
clear Data
Data=Data_sort(:,1:i)
[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,Data,Train_label); %4改好的svm但是有nan  3 libsvm
Result_Data_sort(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc];
end



end
%% index选择
Feature_Show_Num=sum(All_index_feature_selection);
[Feature_Show_Num_sort,filp_index]=sort(Feature_Show_Num,'descend');
Final_index=(1:select_feature_number);
save('Data_featureselection','All_Train_feature_scale','All_Train_label','All_Test_feature_scale','All_Test_label','All_P','All_Order','Feature_Show_Num_sort');

for k=1:k_bootstrap
%% 1/10的数据    测试集
Train_final_featureset_scale=All_Train_feature_scale(:,Final_index,k);
Test_final_featureset_scale=All_Test_feature_scale(:,Final_index,k);
Train_label;
Test_label;

[auc,acc,sens,spec,ppv,npv,mcc]=classify(iteration,flag,Train_final_featureset_scale,Train_label,Test_final_featureset_scale,Test_label );
Result(k,:)=[auc,acc,sens,spec,ppv,npv,mcc]

end
% for i=select_feature_number
% Data=Data_sort(:,1:i);
% [auc,acc,sens,spec,ppv,npv,mcc]=classify(iteration,flag,Data,Train_label,Test_featureselect_scale,Test_label );
% Result(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% 
% end 

Result_final_mean=mean(Result)
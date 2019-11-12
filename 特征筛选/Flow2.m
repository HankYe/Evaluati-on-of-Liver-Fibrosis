% clear all
clc
%十折
%% load data
% load('TI_Features_324p2GAmrmr_.mat')
% load('TI_Label_169.mat')
% load('Feature_TI_System_521_0628.mat')
% load('Label_TI_System_521.mat')
load('yesornoforsample_521to455.mat')
% load('yesornoforsample_521to181')
% load('yesornoforsample_see521to299.mat')
% load('yesornoforsample_521to430.mat')
% load('YesOrNoForPatients.mat') %350
% load('yesornoforsample_506to181.mat')
% load('yesornoforall_1142.mat')
load('Feature_TI_System_521_0626.mat', 'Feature_TI_System_521_0626')
load('Feature_TI_System_521_0626.mat', 'Label_TI_System_521')
% load('matlabdata20160710mat.mat', 'index_313to637_of181')
% load('matlabdata20160710mat.mat', 'index_313to637to1142_of181')
load('date_order_of_in_hospital_521.mat')
load('yesornoforfeature.mat') %1142to637
load('TypeandNamelist_1142.mat')

%% set
%按入院时间排序
All_feature_time=Feature_TI_System_521_0626(date_order_of_in_hospital_521,:);
All_label_time=Label_TI_System_521(date_order_of_in_hospital_521,:);  %Label_TI_System_521为0,1

% yesornoforsample_521to430_2=yesornoforsample_521to430(date_order_of_in_hospital_521,:);
yesornoforsample_521to455_2=yesornoforsample_521to455(date_order_of_in_hospital_521,:);
% yesornoforsample_see521to299_2=yesornoforsample_see521to299(date_order_of_in_hospital_521,:);

All_feature_time2=All_feature_time( find(yesornoforsample_521to455_2==0),find(yesornoforfeature==1)); %521to455 1142to637
All_label=All_label_time(find(yesornoforsample_521to455_2==0),:)+1; %0,1改成1,2

%删除0和缺失值 
index_nan=find(~isnan(mean(All_feature_time2)));
index_zero=find(mean(All_feature_time2)~=0);

index_nanandzero=intersect(index_nan,index_zero) ;
features_nozeroandnan=All_feature_time2(:,index_nanandzero);
All_feature_time3=features_nozeroandnan; 

%637to614   
index_origin=find(yesornoforfeature==1);%637of1142
index_nanandzero2=zeros(1,size(All_feature_time3,2));
for i=1:size(All_feature_time3,2)
    index_nanandzero2(i)=index_origin(index_nanandzero(i));%614of1142
end

%% 参数设置
select_feature_number=75;  %选出特征个数
num_trainandtest=300;  %训练测试集样本数   300 350 400    其余为验证集

for l=75 %15:10:85
    select_feature_number=l;
method_p='ttest2'; %choose from  'ttest2'  'ranksum'
% method_featureselection='SRC'; % SRC SRCnew LASSO mrmr Fisher_score reliefF F_score
j_bootstrap=55; %500
k_flod=10;
iteration=1;  %AdaBoost迭代次数 svm不需要
flag=1;   %使用多分类svm，要求数据分两类排列好
p_threshold=0.05;   %单因素分析p值的阈值

 
[All_feature_scale,A0,B0] = scaling(All_feature_time3);

Trainandtest_feature=All_feature_scale(1:num_trainandtest,:);
Trainandtest_label=All_label(1:num_trainandtest,:);
validation_feature=All_feature_scale(num_trainandtest+1:end,:);
validation_label= All_label (num_trainandtest+1:end,:);

for j=38:j_bootstrap
    %% 数据集构造
    % 特征筛选训练集/训练集  9/10
    % 测试集  1/10
    % 验证集
    % clear Train_feature_scale
    % clear Train_label
    % clear Test_feature
    % clear Test_label
    % clear Train_feature_scale
    % clear Test_feature_scale
    % clear output_test
    % clear P
    % clear H2
    % clear index
    % clear featureselection_index
    
    
    %% 十分数据集的产生
    %----------------------------------------------------------------
    % kflod input:
    % Trainandtest_feature;
    % Trainandtest_label;
    %----------------------------------------------------------------
    [sampletotal,imgCol]=size(Trainandtest_feature);
    
    A = find(Trainandtest_label==1); % malignant index
    B = find(Trainandtest_label==2); % benign index
    
    n_malignant = length(A); % malignant
    n_benign = sampletotal- n_malignant; % benign
    
    in_malignant=fix(n_malignant/10); %integer number of malignant
    rn_malignant=mod(n_malignant,10); %residue number of malignant
    in_benign=fix(n_benign/10); %integer number of benign
    rn_benign=mod(n_benign,10); %residue number of benign
    % n_ir_mn=in_malignant+rn_malignant+in_benign+rn_benign;
    
    A= A(randperm(n_malignant)); %打乱索引A的顺序得到索引Ai
    B= B(randperm(n_benign));
    
    index_train=[];%=zeros(,kflod);
    index_test=[];
    label_kflod=[];
    % label_test_kflod =[];
    index_train2=[];%=zeros(,kflod);
    index_test2=[];
    label_kflod2=[];
    % label_test_kflod2=[];
    
    for time=1:k_flod
        if time<k_flod
            ceiling_malignant=in_malignant*(time-1)+1;
            floor_malignant=time*in_malignant;
            %        index_train(1:in_malignant,:,time) = A(ceiling_malignant:floor_malignant)
            index_train(1:in_malignant,:,time) = A(ceiling_malignant:floor_malignant);
            index_test(1:in_malignant,:,time) = A(ceiling_malignant:floor_malignant);
            
            ceiling_benign=in_benign*(time-1)+1;
            floor_benign=time*in_benign;
            index_train((in_malignant+1):(in_malignant+in_benign),:,time) = B(ceiling_benign:floor_benign);
            index_test((in_malignant+1):(in_malignant+in_benign),:,time) = B(ceiling_benign:floor_benign);
            
        else
            ceiling_malignant=(time-1)*in_malignant+1;
            floor_malignant=n_malignant;
            index_train2(1:(in_malignant+rn_malignant),:) = A(ceiling_malignant:floor_malignant);
            index_test2(1:(in_malignant+rn_malignant),:) = A(ceiling_malignant:floor_malignant);
            
            ceiling_benign=(time-1)*in_benign+1;
            floor_benign=n_benign;
            index_train2((in_malignant+rn_malignant+1):(in_malignant+rn_malignant+in_benign+rn_benign),:) = B(ceiling_benign:floor_benign);
            index_test2((in_malignant+rn_malignant+1):(in_malignant+rn_malignant+in_benign+rn_benign),:) = B(ceiling_benign:floor_benign);
        end
        
    end
    
    
    % Rearranged Target vector, 1 for malignant, 2 for benign
    label_kflod=[ones(1,in_malignant), 1+ones(1,in_benign)]'; %一组的标签[malignant benign]
    label_kflod2=[ones(1,in_malignant+rn_malignant), 1+ones(1,in_benign+rn_benign)]'; %最后一组的标签 [malignant benign]
    
    Train_index=[];
    Train_label=[];
    Test_index=[];
    Test_label=[];
    
    %循环十次
    for k=1:k_flod
        switch k
            case {1},
                Train_index = [index_train(:,:,2)', index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,1);
                Test_label = label_kflod;
            case {2},
                Train_index = [index_train(:,:,1)', index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,2);
                Test_label = label_kflod;
            case {3},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,3);
                Test_label = label_kflod;
            case {4},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,4);
                Test_label = label_kflod;
            case {5},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,5);
                Test_label = label_kflod;
            case {6},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,6);
                Test_label = label_kflod;
            case {7},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,8)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,7);
                Test_label = label_kflod;
            case {8},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,9)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,8);
                Test_label = label_kflod;
            case {9},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train2(:,:)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod2',]';
                Test_index = index_test(:,:,9);
                Test_label = label_kflod;
            case {10},
                Train_index = [index_train(:,:,1)', index_train(:,:,2)',index_train(:,:,3)',index_train(:,:,4)',index_train(:,:,5)',index_train(:,:,6)',index_train(:,:,7)',index_train(:,:,8)',index_train(:,:,9)']';
                Train_label = [label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',label_kflod',]';
                Test_index = index_test2(:,:);
                Test_label = label_kflod2;
        end
        
        numgeneral = size(index_train,1);
        numoutstanding = size(index_train2,1);

        % Train_feature,Train_label,Test_feature,Test_label
        Train_index2=[Train_index(find(Train_label==1),:)' Train_index(find(Train_label==2),:)']';   
        Test_index2=[Test_index(find(Test_label==1),:)' Test_index(find(Test_label==2),:)']';
               
        Train_feature_scale=Trainandtest_feature(Train_index2,:);
        Train_label= Trainandtest_label(Train_index2,:);
        Test_feature_scale=Trainandtest_feature(Test_index2,:);;
        Test_label= Trainandtest_label(Test_index2,:);
%         bootstrap_index_data_Train(k,:,j)=Train_index2; bootstrap_index_data_Test(k,:,j)=Test_index2; %保存bootstrap的索引
        %----------------------------------------------------------------
        % kflod output:
        % Train_feature_scale;
        % Train_label;
        % Test_feature_scale;
        % Test_label;
        %----------------------------------------------------------------
        
        %label 1 2
        % [Train_feature_scale,A0,B0] = scaling(Train_feature);
        % Test_feature_scale = scaling(Test_feature,1,A0,B0); %归一化
        
        % %每次bootstrap十折得到的数据都存为一层数据，存一个整体的数据
        % All_Train_feature_scale(:,:,k)=Train_feature_scale;
        % All_Train_label(:,:,k)=Train_label;
        % All_Test_feature_scale(:,:,k)=Test_feature_scale;
        % All_Test_label(:,:,k)=Test_label;
        
        %% 十分法训练+测试
        %% 9/10的数据    特征筛选训练集
        [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,Train_feature_scale,Train_label); %4改好的svm但是有nan  3 libsvm
        Result_origin(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]  %十分之九用原始特征的效果
        
        %% 单因素分析 p值
        
        [output_test,P,H2,index_p]=p_value_function(Train_feature_scale,Train_label,p_threshold,method_p);

        output_test=Train_feature_scale(:,index_p);
        
%         index_o2p=zeros(1,size(index_p,2));
%         for i=1:size(index_p,1)
%             index_o2p(i)=index_nanandzero2(index_p(i));
%         end
          
%         k_p_index(k,:)=index_p;
        
        [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,output_test,Train_label); %4改好的svm但是有nan  3 libsvm
        Result_p(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]  %十分之九用p值筛选的效果
        
        %% 特征选择
        [~,filp_index,Data_sort]=feature_selection(output_test,Train_label,'SRC');
        % method= GA_mRMR    SRCnew  SRC  LASSO  mrmr  Fisher_score  reliefF  F_score
        
        Data_sort=Train_feature_scale(:,index_feature_selection ) ;
        
        %阈值确定实验
        % for i=5:5:50
        % clear Data
        % Data=Data_sort(:,1:i)
        % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,Data,Train_label); %4改好的svm但是有nan  3 libsvm
        % Result_Data_SRC(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc];
        % end
        
        %SRC对应前75个
        clear Data
        Data=Data_sort(:,1:select_feature_number);
        [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH( 5,1,Data,Train_label); %4改好的svm但是有nan  3 libsvm
        Result_SRC(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
        
        k_index_feature_selection(k,:)=index_feature_selection(:,1:select_feature_number);
        
        %% 1/10的数据    测试集  第十折进行验证
        X_train=Train_feature_scale(:,index_feature_selection );
        X_test=Test_feature_scale(:,index_feature_selection );
        T_train=Train_label;
        T_test=Test_label;
        [classes,testLabel,scores]=solo_validation(X_train,X_test,T_train,T_test );
        
        %记下每次十折的预测值
        if k<10
            Label_predict((k-1)*numgeneral+1:k*numgeneral,:) = classes(:,:); % test result
            Label_golden((k-1)*numgeneral+1:k*numgeneral,:)=testLabel(:,:);
            Scores((k-1)*numgeneral+1:k*numgeneral,:)=scores(:,:);
        else
            Label_predict((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:) = classes(:,:); % test result
            Label_golden((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=testLabel(:,:);
            Scores((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=scores(:,:);
        end
    i,
    j,
    end  %kflod's end
    
    %% 第j次bootstrap的十次十折预测结果的评价
    [ auc,acc,sens,spec,ppv,npv,mcc]= assess(Label_predict,Label_golden,Scores);
    Result_test(j,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc] %第j次bootstrap的用选出特征的十分之一的最终结果
    
    %第j次bootstrap的每次特征筛选的结果的平均值
    Result_origin2(j,:)=mean(Result_origin);
    Result_p2(j,:)=mean( Result_p);
    Result_gamrmr_output2(j,:)=mean( Result_gamrmr_output);
    Result_SRC2(j,:)=mean(Result_SRC);
    
    temp_best=zeros(k_flod,imgCol);
    for i=1:k_flod
        temp_best(i,k_index_feature_selection(i,:))=1;
    end
    union=sum(temp_best);
    [~,k_final_index]=sort(union,'descend');
    
    j_index_feature_selection(j,:)=k_final_index(:,1:select_feature_number);
    
end  %bootstrap's end
%训练和测试的平均结果
Result_trainandtest_bootstrap_mean=mean( Result_test)
Result_origin_bootstrap_mean =mean(Result_origin2);
Result_p_bootstrap_mean =mean( Result_p2);
Result_gamrmr_output_bootstrap_mean=mean( Result_gamrmr_output2);
Result_SRC_bootstrap_mean=mean(Result_SRC2);
    
%% index选择
temp_best=zeros(j_bootstrap,imgCol);
for i=1:j_bootstrap
    temp_best(i,j_index_feature_selection(i,:))=1;
end
union=sum(temp_best);
[~,j_final_index]=sort(union,'descend');

Final_index_feature_selection=j_final_index(:,1:select_feature_number);

%% 全新验证集
for i=1:select_feature_number
X_train_final=Trainandtest_feature(:,Final_index_feature_selection(:,1:i));
T_train_final=Trainandtest_label;
X_test_final=validation_feature(:,Final_index_feature_selection(:,1:i) );
T_test_final=validation_label;
[classes,testLabel,scores]=solo_validation(X_train_final,X_test_final,T_train_final,T_test_final );
clear Label_predict
clear Label_golden
clear Scores
Label_predict(:,:) = classes(:,:); % test result
Label_golden(:,:)=testLabel(:,:);
Scores(:,:)=scores(:,:);
[ auc,acc,sens,spec,ppv,npv,mcc]= assess(Label_predict,Label_golden,Scores);
Result_validation(i,:)=[auc,acc,sens,spec,ppv,npv,mcc]  %验证集分类结果
end
[best_AUC_validation,index_best_AUC_position]=max(Result_validation(:,1))
Result_validation_final=Result_validation(index_best_AUC_position,:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%637to614  %%614to1142 
Final_index_feature_selection_1142=zeros(1,size(Final_index_feature_selection,2));
for i=1:size(Final_index_feature_selection,2)
    Final_index_feature_selection_1142(i)=index_nanandzero2(Final_index_feature_selection(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_Result_origin_bootstrap_mean(l,:) =Result_origin_bootstrap_mean;
F_Result_p_bootstrap_mean(l,:) =Result_p_bootstrap_mean ;
F_Result_gamrmr_output_bootstrap_mean(l,:)=Result_gamrmr_output_bootstrap_mean;
F_Result_SRC_bootstrap_mean(l,:)=Result_SRC_bootstrap_mean;

F_Result_trainandtest_bootstrap_mean(l,:)=Result_trainandtest_bootstrap_mean;
F_Result_validation_final(l,:)=Result_validation_final;
save('data_20170728','F_Result_trainandtest_bootstrap_mean','F_Result_validation_final','F_Result_origin_bootstrap_mean','F_Result_p_bootstrap_mean','F_Result_gamrmr_output_bootstrap_mean','F_Result_SRC_bootstrap_mean');
end

% save('','Final_index_feature_selection','Result_validation_final','Result_validation','','','','','','','','','','','','','','');
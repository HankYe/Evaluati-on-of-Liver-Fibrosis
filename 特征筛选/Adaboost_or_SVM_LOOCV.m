function [auc,acc,sens,spec,ppv,npv,mcc,Scores]=Adaboost_or_SVM_LOOCV(iteration,flag,X_scale,T_scale)
% global L;
%  ADABOOST_and_SVM_train_and_test输入参数解释
%  C,gamma,iteration,flag,X_train,X_test,T_scale
%  C:SVM的错误代价系数，默认为1
%  gamma:SVM的径向基核函数的参数 gamma，默认值为1/输入特征总数
%  iteration为Adaboost的的迭代次数
%  flag为选择分类器是SVM还是Adaboost
%  X_train为训练集，同X_test
%  因为此函数用LOOCV留一法，输入时训练集和测试集一样，但是实际上是留一法，
%  每次取出来一个做测试，其他的都是训练集。之后再循环，直到每个都做过一次测试集
%  T_scale为数据的label，尤0和1组成。但是在程序运行时会被转化为1,2的label
%  需要设置空文件夹 程序所在文件夹目录\MAT\  用来存放中间数据

% 整个训练过程，包括读入图像，多类 SVM 训练或者多类 Adaboost，各个阶段的处理结果分别保存至文件：
%   将 scaling 的各维上、下界信息保存至 Mat/scaling.mat
%   将 PCA 降维并且 scaling 后的数据保存至 Mat/trainData.mat
%   将多类 SVM 的训练信息保存至 Mat/multiSVMTrain.mat
%   将多类 Adaboost 的训练信息保存至 Mat/multiAdaboostTrain.mat

%  此程序由Adaboost人脸识别程序改编而来，删除了数据读入、PCA降维部分

%demo：ADABOOST_and_SVM_train_and_test(1,1/199,8,1,Features,Features,IDH1); %使用SVM，特征为88*199


X_train=X_scale;
X_test=X_scale;

% display(' ');
% display(' ');
% display('训练开始...');

%% 数据载入和转换格式
% C=C_sel;
% gamma=g_sel;

[imgRow,imgCol]=size(X_train);
C=1;
gamma=1/imgCol;
if (flag==5)
    for i=1:imgRow   %0,1的label换为1，2，3...形式
        if(T_scale(i)==1)
            T_scale(i)=0;
        else T_scale(i)=1;
        end
    end
else
    T_scale=T_scale+1;
%     for i=1:imgRow   %0,1的label换为1，2，3...形式
%         if(T_scale(i)==1)
%             T_scale(i)=2;
%         else T_scale(i)=1;
%         end
%     end
end
nLabel=2;%分类数目

T=T_scale;
%% 选出label对应的例子，分前后两部分
e=1e-10;
sampletotal=size(X_train,1);
if (flag==5)
    M = find(T==0); % malignant index 1-88 
    B = find(T==1); % benign index 1-88
else
    M = find(T==1); % malignant index 1-88 
    B = find(T==2); % benign index 1-88
end
n_malignant = length(M); % malignant number
n_benign = sampletotal- n_malignant; % benign number

% Rearrange X M和B的拼起来
% Mi= M(randperm(n_malignant)); %打乱索引M的顺序得到索引Mi 1-88的63个
X_train_scale(1:n_malignant,:) = X_train(M,:);
X_test_scale(1:n_malignant,:) = X_test(M,:);

% Bi= B(randperm(n_benign));
X_train_scale(n_malignant+1:sampletotal,:) = X_train(B,:);
X_test_scale(n_malignant+1:sampletotal,:) = X_test(B,:);

if (flag==5)
    T_scale=[1-ones(1,n_malignant), ones(1,n_benign)]'; % [malignant benign]
else
    % Rearranged Target vector, 1 for malignant, 2 for benign
    T_scale=[ones(1,n_malignant), 1+ones(1,n_benign)]'; % [malignant benign]
end
% X_train_scale = X_train;
% X_test_scale = X_test;
% T_scale=T;

%% LOOCV的 ADABOOST or SVM 程序主体

for k = 1:sampletotal
%% 留一法构造训练集测试集   
    if k == 1
        D_train = X_train_scale(2:sampletotal,:);
        T_train = T_scale(2:sampletotal);
%         %%%
%                 D_test = X_test_scale(1:3,:);
%                 T_test = T_scale(1:3);

    elseif k == sampletotal
        D_train = X_train_scale(1:sampletotal-1,:);
        T_train = T_scale(1:sampletotal-1);    
       
%         %%%
%                 %%%
%                 D_test = X_test_scale(sampletotal-3:sampletotal,:);
%                 T_test = T_scale(sampletotal-3:sampletotal);
    else
        D_train = [X_train_scale(1:k-1,:)' X_train_scale(k+1:sampletotal,:)']';
        T_train = [T_scale(1:k-1)' T_scale(k+1:sampletotal)']';
        
%         %%%%%
%         D_test = X_test_scale(k:k+3,:);
%         T_test = T_scale(k:k+3);
    end
    
    D_test = X_test_scale(k,:);
    T_test = T_scale(k);

    
%% 训练部分 
X = D_train;
dataLabel= T_train;
nLabel=2;

% display('归一化开始...');
% display('.........');
[X,A0,B0] = scaling(X);
% save('Mat/scaling.mat', 'A0', 'B0');
% 保存 scaling 后的训练数据至 trainData.mat

TrainData = X;
trainLabel = dataLabel;
% save('Mat/trainData.mat', 'TrainData', 'trainLabel');
% display('归一化完成...');

% display('训练分类器开始，这个过程可能会花上几分钟.........................');
if(flag==1) %SVM模式 
    for iLabel = 1:nLabel
        nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %存储不同label例子的个数
    end

    multiSVMStruct = multiSVMTrain(TrainData, nSplPerClass, nLabel, C, gamma);
    %参数意义：训练集 label各种类的例子个数 label值 c gamma
%     display('正在保存SVM训练结果...');
%     save('Mat/multiSVMTrain.mat', 'multiSVMStruct');
else if(flag==0) %adaboost模式 
    for iLabel = 1:nLabel
          nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %存储不同label例子的个数
    end
    
    multiAdaboostStruct = multiAdaboostTrain(TrainData, nSplPerClass ,nLabel,iteration);
    %参数意义：训练集 label各种类的例子个数 label值 迭代次数（根据不同数据选择不同迭代次数）
%     display('正在保存Adaboost训练结果...');
%     save('Mat/multiAdaboostTrain.mat', 'multiAdaboostStruct');
elseif (flag==5)
        a =glmfit(TrainData,trainLabel,'binomial', 'link', 'logit');
        
    else
        display('flag 错误，请输入0 or 1，对应svm和adaboost');
    end
end
% display('..............................');
% display('训练完成。');

%% 测试部分
%     display(' ');
%     display(' ');
%     display('测试开始...');
%     display('载入训练参数...');
     
 if(flag==1)
    %数据准备
    TestData=  D_test;
    testLabel= T_test;
    TestData = scaling(TestData,1,A0,B0); %归一化
%     display('..............................');

    % 多类 SVM 分类
%     display('测试集识别中...');
    [classes,scores] = multiSVMClassify(TestData, multiSVMStruct);
%     display('..............................');
    
    % 计算识别率
    nError = sum(classes ~= testLabel);
    acc = 1 - nError/length(testLabel);
%     display(['LOOCV+SVM对于测试集的识别率为', num2str(acc*100), '%','第',num2str(k),'次']);
%     Accuracy(k)=accuracy;
   
    Label_predict(k) = classes; % test result
    Label_golden(k)=testLabel';
    Scores(k)=scores;
    
 else if(flag==0)
    %数据准备
    TestData= D_test;
    testLabel= T_test;
    TestData = scaling(TestData,1,A0,B0);
%     display('..............................');
    
     % 多类 AdaBoost 分类
%     display('测试集识别中...');
    [classes,scores] = multiAdaboostClassify(TestData,multiAdaboostStruct);
    display('..............................');
    
    % 计算识别率
    nError = sum(classes ~= testLabel);
    acc = 1 - nError/length(testLabel);
%     display(['LOOCV+Adaboost对于测试集的识别率为', num2str(acc*100), '%' '第',num2str(k),'次']);
%     Accuracy(k)=accuracy;
    Label_predict(k) = classes; % test result
    Label_golden(k)=testLabel';   
    Scores(k)=scores(:,1);
%     Scores(k)=L;
     elseif (flag==5)
         %数据准备
         TestData= D_test;
         testLabel= T_test;
         TestData = scaling(TestData,1,A0,B0);
%          display('..............................');
         
         % logistics
%          display('测试集识别中...');
         
         [logitFit] = glmval(a,TestData,'logit');
         scores=logitFit;
         classes=sign(logitFit-0.5)+1;
%          display('..............................');
         
         % 计算识别率
         nError = sum(classes ~= testLabel);
         acc = 1 - nError/length(testLabel);
%          display(['LOOCV+logistics对于测试集的识别率为', num2str(acc*100), '%' '第',num2str(k),'次']);
         %     Accuracy(k)=accuracy;
         Label_predict(k) = classes; % test result
         Label_golden(k)=testLabel';
         Scores(k)=scores(:,1);
         %     Scores(k)=L;
   
     else
         display('flag 错误，请输入0或1 ');
     end
end
% Accuracy_singel(i)=mean(Accuracy)
end
 
%% 判断是否分错了 分错了就显示1，分对了为0。用于研究原始数据是否准确，反馈医生
for i=1:imgRow
    if( Label_predict(i)~=Label_golden(i))
        TEMP(i)=1;
    else
        TEMP(i)=0;
    end
    
end

%% 评价指标
[ auc,acc,sens,spec,ppv,npv,mcc]= assess(Label_predict,Label_golden,Scores);
 
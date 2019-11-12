function [auc,acc,sens,spec,ppv,npv,mcc,Scores]=Adaboost_or_SVM_LOOCV(iteration,flag,X_scale,T_scale)
% global L;
%  ADABOOST_and_SVM_train_and_test�����������
%  C,gamma,iteration,flag,X_train,X_test,T_scale
%  C:SVM�Ĵ������ϵ����Ĭ��Ϊ1
%  gamma:SVM�ľ�����˺����Ĳ��� gamma��Ĭ��ֵΪ1/������������
%  iterationΪAdaboost�ĵĵ�������
%  flagΪѡ���������SVM����Adaboost
%  X_trainΪѵ������ͬX_test
%  ��Ϊ�˺�����LOOCV��һ��������ʱѵ�����Ͳ��Լ�һ��������ʵ��������һ����
%  ÿ��ȡ����һ�������ԣ������Ķ���ѵ������֮����ѭ����ֱ��ÿ��������һ�β��Լ�
%  T_scaleΪ���ݵ�label����0��1��ɡ������ڳ�������ʱ�ᱻת��Ϊ1,2��label
%  ��Ҫ���ÿ��ļ��� ���������ļ���Ŀ¼\MAT\  ��������м�����

% ����ѵ�����̣���������ͼ�񣬶��� SVM ѵ�����߶��� Adaboost�������׶εĴ������ֱ𱣴����ļ���
%   �� scaling �ĸ�ά�ϡ��½���Ϣ������ Mat/scaling.mat
%   �� PCA ��ά���� scaling ������ݱ����� Mat/trainData.mat
%   ������ SVM ��ѵ����Ϣ������ Mat/multiSVMTrain.mat
%   ������ Adaboost ��ѵ����Ϣ������ Mat/multiAdaboostTrain.mat

%  �˳�����Adaboost����ʶ�����ı������ɾ�������ݶ��롢PCA��ά����

%demo��ADABOOST_and_SVM_train_and_test(1,1/199,8,1,Features,Features,IDH1); %ʹ��SVM������Ϊ88*199


X_train=X_scale;
X_test=X_scale;

% display(' ');
% display(' ');
% display('ѵ����ʼ...');

%% ���������ת����ʽ
% C=C_sel;
% gamma=g_sel;

[imgRow,imgCol]=size(X_train);
C=1;
gamma=1/imgCol;
if (flag==5)
    for i=1:imgRow   %0,1��label��Ϊ1��2��3...��ʽ
        if(T_scale(i)==1)
            T_scale(i)=0;
        else T_scale(i)=1;
        end
    end
else
    T_scale=T_scale+1;
%     for i=1:imgRow   %0,1��label��Ϊ1��2��3...��ʽ
%         if(T_scale(i)==1)
%             T_scale(i)=2;
%         else T_scale(i)=1;
%         end
%     end
end
nLabel=2;%������Ŀ

T=T_scale;
%% ѡ��label��Ӧ�����ӣ���ǰ��������
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

% Rearrange X M��B��ƴ����
% Mi= M(randperm(n_malignant)); %��������M��˳��õ�����Mi 1-88��63��
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

%% LOOCV�� ADABOOST or SVM ��������

for k = 1:sampletotal
%% ��һ������ѵ�������Լ�   
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

    
%% ѵ������ 
X = D_train;
dataLabel= T_train;
nLabel=2;

% display('��һ����ʼ...');
% display('.........');
[X,A0,B0] = scaling(X);
% save('Mat/scaling.mat', 'A0', 'B0');
% ���� scaling ���ѵ�������� trainData.mat

TrainData = X;
trainLabel = dataLabel;
% save('Mat/trainData.mat', 'TrainData', 'trainLabel');
% display('��һ�����...');

% display('ѵ����������ʼ��������̿��ܻỨ�ϼ�����.........................');
if(flag==1) %SVMģʽ 
    for iLabel = 1:nLabel
        nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %�洢��ͬlabel���ӵĸ���
    end

    multiSVMStruct = multiSVMTrain(TrainData, nSplPerClass, nLabel, C, gamma);
    %�������壺ѵ���� label����������Ӹ��� labelֵ c gamma
%     display('���ڱ���SVMѵ�����...');
%     save('Mat/multiSVMTrain.mat', 'multiSVMStruct');
else if(flag==0) %adaboostģʽ 
    for iLabel = 1:nLabel
          nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %�洢��ͬlabel���ӵĸ���
    end
    
    multiAdaboostStruct = multiAdaboostTrain(TrainData, nSplPerClass ,nLabel,iteration);
    %�������壺ѵ���� label����������Ӹ��� labelֵ �������������ݲ�ͬ����ѡ��ͬ����������
%     display('���ڱ���Adaboostѵ�����...');
%     save('Mat/multiAdaboostTrain.mat', 'multiAdaboostStruct');
elseif (flag==5)
        a =glmfit(TrainData,trainLabel,'binomial', 'link', 'logit');
        
    else
        display('flag ����������0 or 1����Ӧsvm��adaboost');
    end
end
% display('..............................');
% display('ѵ����ɡ�');

%% ���Բ���
%     display(' ');
%     display(' ');
%     display('���Կ�ʼ...');
%     display('����ѵ������...');
     
 if(flag==1)
    %����׼��
    TestData=  D_test;
    testLabel= T_test;
    TestData = scaling(TestData,1,A0,B0); %��һ��
%     display('..............................');

    % ���� SVM ����
%     display('���Լ�ʶ����...');
    [classes,scores] = multiSVMClassify(TestData, multiSVMStruct);
%     display('..............................');
    
    % ����ʶ����
    nError = sum(classes ~= testLabel);
    acc = 1 - nError/length(testLabel);
%     display(['LOOCV+SVM���ڲ��Լ���ʶ����Ϊ', num2str(acc*100), '%','��',num2str(k),'��']);
%     Accuracy(k)=accuracy;
   
    Label_predict(k) = classes; % test result
    Label_golden(k)=testLabel';
    Scores(k)=scores;
    
 else if(flag==0)
    %����׼��
    TestData= D_test;
    testLabel= T_test;
    TestData = scaling(TestData,1,A0,B0);
%     display('..............................');
    
     % ���� AdaBoost ����
%     display('���Լ�ʶ����...');
    [classes,scores] = multiAdaboostClassify(TestData,multiAdaboostStruct);
    display('..............................');
    
    % ����ʶ����
    nError = sum(classes ~= testLabel);
    acc = 1 - nError/length(testLabel);
%     display(['LOOCV+Adaboost���ڲ��Լ���ʶ����Ϊ', num2str(acc*100), '%' '��',num2str(k),'��']);
%     Accuracy(k)=accuracy;
    Label_predict(k) = classes; % test result
    Label_golden(k)=testLabel';   
    Scores(k)=scores(:,1);
%     Scores(k)=L;
     elseif (flag==5)
         %����׼��
         TestData= D_test;
         testLabel= T_test;
         TestData = scaling(TestData,1,A0,B0);
%          display('..............................');
         
         % logistics
%          display('���Լ�ʶ����...');
         
         [logitFit] = glmval(a,TestData,'logit');
         scores=logitFit;
         classes=sign(logitFit-0.5)+1;
%          display('..............................');
         
         % ����ʶ����
         nError = sum(classes ~= testLabel);
         acc = 1 - nError/length(testLabel);
%          display(['LOOCV+logistics���ڲ��Լ���ʶ����Ϊ', num2str(acc*100), '%' '��',num2str(k),'��']);
         %     Accuracy(k)=accuracy;
         Label_predict(k) = classes; % test result
         Label_golden(k)=testLabel';
         Scores(k)=scores(:,1);
         %     Scores(k)=L;
   
     else
         display('flag ����������0��1 ');
     end
end
% Accuracy_singel(i)=mean(Accuracy)
end
 
%% �ж��Ƿ�ִ��� �ִ��˾���ʾ1���ֶ���Ϊ0�������о�ԭʼ�����Ƿ�׼ȷ������ҽ��
for i=1:imgRow
    if( Label_predict(i)~=Label_golden(i))
        TEMP(i)=1;
    else
        TEMP(i)=0;
    end
    
end

%% ����ָ��
[ auc,acc,sens,spec,ppv,npv,mcc]= assess(Label_predict,Label_golden,Scores);
 
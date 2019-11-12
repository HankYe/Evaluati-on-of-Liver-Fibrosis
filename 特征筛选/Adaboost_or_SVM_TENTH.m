% function [auc,acc,sens,spec,ppv,npv,mcc,Scores]=Adaboost_or_SVM_TENTH(iteration,flag,X_scale,T_scale)

function [auc,acc,sens,spec,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(iteration,flag,X_scale,T_scale)
% display(' ');
% display(' ');
% display('ѵ����ʼ...');

%% ���������ת����ʽ
% C=C_sel;
% gamma=g_sel;
[imgRow,imgCol]=size(X_scale);
C=1;
gamma=1/imgCol;

for i=1:imgRow   %0,1��label��Ϊ1��2��3...��ʽ
if(T_scale(i)==1)
    T_scale(i)=2;
else T_scale(i)=1;
end
end

nLabel=2;%������Ŀ

T=T_scale;
%% ѡ��label��Ӧ�����ӣ���ǰ��������
e=1e-10;
sampletotal=size(X_scale,1);

M = find(T==1); % malignant index 
B = find(T==2); % benign index 

n_malignant = length(M); % malignant
n_benign = sampletotal- n_malignant; % benign

in_malignant=fix(n_malignant/10); %integer number of malignant
rn_malignant=mod(n_malignant,10); %residue number of malignant
in_benign=fix(n_benign/10); %integer number of benign
rn_benign=mod(n_benign,10); %residue number of benign
n_ir_mn=in_malignant+rn_malignant+in_benign+rn_benign;

% M= M(randperm(n_malignant)); %��������M��˳��õ�����Mi 1-88��63��
% B= B(randperm(n_benign));

for time=1:10
    if time<10
        % Rearrange X M��B��ƴ����
        ceiling_malignant=in_malignant*(time-1)+1;
        floor_malignant=time*in_malignant;
        X_train_scale(1:in_malignant,:,time) = X_scale(M(ceiling_malignant:floor_malignant),:);
        X_test_scale(1:in_malignant,:,time) = X_scale(M(ceiling_malignant:floor_malignant),:);
    
        ceiling_benign=in_benign*(time-1)+1;
        floor_benign=time*in_benign;
        X_train_scale((in_malignant+1):(in_malignant+in_benign),:,time) = X_scale(B(ceiling_benign:floor_benign),:);
        X_test_scale((in_malignant+1):(in_malignant+in_benign),:,time) = X_scale(B(ceiling_benign:floor_benign),:);
   else 
        ceiling_malignant=(time-1)*in_malignant+1;
        floor_malignant=n_malignant;
        X_train_scale2(1:(in_malignant+rn_malignant),:) = X_scale(M(ceiling_malignant:floor_malignant),:);
        X_test_scale2(1:(in_malignant+rn_malignant),:) = X_scale(M(ceiling_malignant:floor_malignant),:);
   
        ceiling_benign=(time-1)*in_benign+1;
        floor_benign=n_benign;
        X_train_scale2((in_malignant+rn_malignant+1):(in_malignant+rn_malignant+in_benign+rn_benign),:) = X_scale(B(ceiling_benign:floor_benign),:);
        X_test_scale2((in_malignant+rn_malignant+1):(in_malignant+rn_malignant+in_benign+rn_benign),:) = X_scale(B(ceiling_benign:floor_benign),:);

    end
end

% Rearranged Target vector, 1 for malignant, 2 for benign
T_scale=[ones(1,in_malignant), 1+ones(1,in_benign)]'; % [malignant benign]
T_scale2=[ones(1,in_malignant+rn_malignant), 1+ones(1,in_benign+rn_benign)]'; % [malignant benign]

D_train=[];
T_train=[];

for k = 1:10
    switch k
        case {1},
            D_train = [X_train_scale(:,:,2)', X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,1);
            T_test = T_scale;
        case {2},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,2);
            T_test = T_scale;
        case {3},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,3);
            T_test = T_scale;
        case {4},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,4);
            T_test = T_scale;
        case {5},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,5);
            T_test = T_scale;            
        case {6},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,6);
            T_test = T_scale;               
        case {7},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,8)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,7);
            T_test = T_scale;       
        case {8},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,9)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,8);
            T_test = T_scale;              
        case {9},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale2(:,:)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale2',]';     
            D_test = X_test_scale(:,:,9);
            T_test = T_scale;              
        case {10},
            D_train = [X_train_scale(:,:,1)', X_train_scale(:,:,2)',X_train_scale(:,:,3)',X_train_scale(:,:,4)',X_train_scale(:,:,5)',X_train_scale(:,:,6)',X_train_scale(:,:,7)',X_train_scale(:,:,8)',X_train_scale(:,:,9)']';
            T_train = [T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',T_scale',]';     
            D_test = X_test_scale2(:,:);
            T_test = T_scale2;             
    end  


D_train2=[D_train(find(T_train==1),:)' D_train(find(T_train==2),:)']';
L_N=length(  find(T_train==1) ) ;
L_P=length(  find(T_train==2) );
T_train2=[ ones(1, L_N  )  1+ones(1,L_P)  ]';  

D_test2=[D_test(find(T_test==1),:)' D_test(find(T_test==2),:)']';
L_N=length(  find(T_test==1) ) ;
L_P=length(  find(T_test==2) );
T_test2=[ ones(1, L_N  )  1+ones(1,L_P)  ]';  

D_train=D_train2;
T_train=T_train2;
D_test=D_test2;
T_test=T_test2;

numgeneral = size(X_train_scale,1);
numoutstanding = size(X_train_scale2,1);
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
if(flag==1) %�����SVMģʽ��Ҫ������ǰһ���һ��ķ���ţ����ܴ����
    for iLabel = 1:nLabel
        nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %�洢��ͬlabel���ӵĸ���
    end

    multiSVMStruct = multiSVMTrain(TrainData, nSplPerClass, nLabel, C, gamma);
    %�������壺ѵ���� label����������Ӹ��� labelֵ c gamma
%     display('���ڱ���SVMѵ�����...');
%     save('Mat/multiSVMTrain.mat', 'multiSVMStruct');
    
    
elseif (flag==0) %adaboostģʽ 
    for iLabel = 1:nLabel
          nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %�洢��ͬlabel���ӵĸ���
    end
    
    multiAdaboostStruct = multiAdaboostTrain(TrainData, nSplPerClass ,nLabel,iteration);
    %�������壺ѵ���� label����������Ӹ��� labelֵ �������������ݲ�ͬ����ѡ��ͬ����������
%     display('���ڱ���Adaboostѵ�����...');
%     save('Mat/multiAdaboostTrain.mat', 'multiAdaboostStruct');
%     else
%         display('flag ����������0 or 1����Ӧsvm��adaboost');
%     end
elseif (flag==3) %libsvm
%             for iLabel = 1:nLabel
%         nSplPerClass(iLabel) = sum( (trainLabel == iLabel) ); %�洢��ͬlabel���ӵĸ���
%             end

            svmset =[ ' -c ',num2str(C), ' -g ', num2str(gamma), ' -b 1'];
%               svmset =[ '-t 2', ' -c ',num2str(C), ' -g ', num2str(gamma), ' -b 1'];
%         svmset = ['-c ',num2str(C),' -g ',num2str(gamma),' -t 3',' -b 1'] ;
%         svmset = ['-c ',num2str(C),' -t 3',' -b 1'] ;
          %�������壺ѵ���� label����������Ӹ��� labelֵ c gamma
        LIBSVMStruct = svmtrain(trainLabel,TrainData, svmset);
%         display('���ڱ���LIBSVMѵ�����...');
%      %�������壺ѵ���� label����������Ӹ��� labelֵ c gamma
%     save('Mat/LIBSVMStruct.mat', 'LIBSVMStruct');
%     end
elseif (flag==4)
    SVMStruct= svm_train(TrainData, trainLabel, 'Kernel_Function', @(TrainData,trainLabel) kfun_rbf(TrainData,trainLabel,gamma), 'boxconstraint', C );
%       SVMStruct= svmtrain(TrainData, trainLabel, 'Kernel_Function', @(TrainData,trainLabel) kfun_rbf(TrainData,trainLabel,gamma));
% % %       SVMStruct= svmtrain(TrainData, trainLabel,'showplot',true);

      %       svmStruct = svmtrain(X(P.training,:),Y(P.training),'showplot',true);
%       C = svmclassify(svmStruct,X(P.test,:),'showplot',true);

else
       display('flag ����������0 or 1 3����Ӧsvm��adaboost  liabsvm');
end
          
%          
% end

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

%     % ����ʶ����
%     nError = sum(classes ~= testLabel);
%     accuracy = 1 - nError/length(testLabel);
% %     display(['ʮ�ַ�TENTH+SVM���ڲ��Լ���ʶ����Ϊ', num2str(accuracy*100), '%','��',num2str(k),'��']);
%     Accuracy(k)=accuracy;
    
    if k<10
        Label_predict((k-1)*numgeneral+1:k*numgeneral,:) = classes(:,:); % test result
        Label_golden((k-1)*numgeneral+1:k*numgeneral,:)=testLabel(:,:);
        Scores((k-1)*numgeneral+1:k*numgeneral,:)=scores(:,:);
    else
        Label_predict((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:) = classes(:,:); % test result
        Label_golden((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=testLabel(:,:);
        Scores((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=scores(:,:);
    end
    
 elseif (flag==3)
    %����׼��
    TestData=  D_test;
    testLabel= T_test;
    TestData = scaling(TestData,1,A0,B0); %��һ��
%     display('..............................');

    % libSVM ����
%     display('���Լ�ʶ����...');
%     display('..............................');
%      
%        display('���Լ�ʶ����...');
       [classes, accuracy, scores] = svmpredict(T_test, D_test,LIBSVMStruct,'-b 1');
%        [classes, accuracy, scores] = svmpredict(testLabel, TestData,LIBSVMStruct,'-b 1');
%      [classes, accuracy, scores] = svmpredict(testLabel, TestData,LIBSVMStruct,'-b 1');
     %    [classes,scores] = multiSVMClassify(TestData, multiSVMStruct);
%      display('..............................');
     
%     % ����ʶ����
%     nError = sum(classes ~= testLabel);
%     accuracy = 1 - nError/length(testLabel);
% %     display(['ʮ�ַ�TENTH+libSVM���ڲ��Լ���ʶ����Ϊ', num2str(accuracy*100), '%','��',num2str(k),'��']);
%     Accuracy(k)=accuracy;
    
    if k<10
        Label_predict((k-1)*numgeneral+1:k*numgeneral,:) = classes(:,:); % test result
        Label_golden((k-1)*numgeneral+1:k*numgeneral,:)=testLabel(:,:);
        Scores((k-1)*numgeneral+1:k*numgeneral,:)=scores(:,:);
    else
        Label_predict((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:) = classes(:,:); % test result
        Label_golden((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=testLabel(:,:);
        Scores((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=scores(:,:);
    end
    
%      end
 elseif (flag==4)
         TestData=  D_test;
         testLabel= T_test;
         TestData = scaling(TestData,1,A0,B0); %��һ��
%          display('svm..............................');
                  
         [classes,scores]  = svmclassify(SVMStruct,TestData);
         
         if k<10
             Label_predict((k-1)*numgeneral+1:k*numgeneral,:) = classes(:,:); % test result
             Label_golden((k-1)*numgeneral+1:k*numgeneral,:)=testLabel(:,:);
             Scores((k-1)*numgeneral+1:k*numgeneral,:)=scores(:,:);
         else
             Label_predict((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:) = classes(:,:); % test result
             Label_golden((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=testLabel(:,:);
             Scores((k-1)*numgeneral+1:(k-1)*numgeneral+numoutstanding,:)=scores(:,:);
         end
        
%              nError = sum(classes ~= testLabel);
%     accuracy = 1 - nError/length(testLabel);
% %     display(['ʮ�ַ�TENTH+svm���ڲ��Լ���ʶ����Ϊ', num2str(accuracy*100), '%' '��',num2str(k),'��']);
%     Accuracy(k)=accuracy;
    
         
 elseif(flag==0)
    %����׼��
    TestData=  D_test;
    testLabel=  T_test;
    TestData = scaling(TestData,1,A0,B0);
%     display('..............................');
    
     % ���� AdaBoost ����
%     display('���Լ�ʶ����...');
    classes = multiAdaboostClassify(TestData,multiAdaboostStruct);
%     display('..............................');
    
%     % ����ʶ����
%     nError = sum(classes ~= testLabel);
%     accuracy = 1 - nError/length(testLabel);
%     display(['ʮ�ַ�TENTH+Adaboost���ڲ��Լ���ʶ����Ϊ', num2str(accuracy*100), '%' '��',num2str(k),'��']);
%     Accuracy(k)=accuracy;
    
    
    
else
        display('flag ����������0��1 ');
end

end
 

% Acc=mean(Accuracy)
% %% �ж��Ƿ�ִ��� �ִ��˾���ʾ1���ֶ���Ϊ0�������о�ԭʼ�����Ƿ�׼ȷ������ҽ��
% for i=1:imgRow
%     if( Label_predict(i)~=Label_golden(i))
%         TEMP(i)=1;
%     else
%         TEMP(i)=0;
%     end
%     
% end
% 
% %% ����ָ��
[ auc,acc,sens,spec,ppv,npv,mcc]= assess(Label_predict,Label_golden,Scores);
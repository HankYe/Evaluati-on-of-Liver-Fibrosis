% %TTEST
% % [h,p,ci,stats] = ttest2(Feature_BC_label1,Feature_BC_label2,0.05)

% load('Feature_TI_System_N12345_953_0531.mat')
% load('output_R2GA_selesct_TI953_N12345.mat')
% load('Feature_TI_System_521_0626.mat')
% load('Feature_TI_System_521_0628.mat')
% load('Label_TI_System_521.mat')
% %特征重构
clear output_test
clear input_b
clear input_m
clear input
clear P
clear H
clear H2
% % load label_1230.mat
% % load feature_1230_4.mat
% load('TI_Label_169.mat')
% load('TI_Feature_SystemV1_950_169.mat')
% 
% I=[1:13 17:19];
% clear feature_bytype; 
% feature_bytype=[];
% for k=I
% %     features=[Feature_TI_System_521_0626(:,1:207)';Feature_TI_System_521_0626(:,395:end)']';
%     features= features_nozeroandnan;
% % % %     feature_bytype(:,:)=features(:,find(Type_nozeroandnan(2,:)==k));
% %         feature_bytype(:,:)=features(:,k);
% 
%     feature_bytype=[feature_bytype features(:,find(Type_nozeroandnan(2,:)==k)) ];
%     
%     
% end
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,feature_bytype,label);
% ResultType(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]; 

% features=Feature_TI_System_521_0626_955(find(YesOrNoForPatients==1),:);%output_test;%feature;%gamrmr_output;%Feature_TI_System_N12345_953_0531
% label=Label_TI_System_521(find(YesOrNoForPatients==1),:);
features=TI_Feature_sift_521;;%(1:169,:);
label=Label_TI_System_521;%(1:169,:);
I=mean(features);
index_NoZeroFeature=find(I~=0);
index_NoNanFeatures=find(~isnan(I));
index_NoZeroAndNanFeatures=intersect(index_NoZeroFeature ,index_NoNanFeatures);
features_nozeroandnan=features(:,index_NoZeroAndNanFeatures);
% Type_nozeroandnan=Type_955(:,index_NoZeroAndNanFeatures);
% TypeandName_nozeroandnan=TypeandName_955(:,index_NoZeroAndNanFeatures);




input=features_nozeroandnan;
% label=Label_TI_System_521(1:169,:);

i_b=1;
i_m=1;

[n_x,n_y]=size(input);

% Feature normalization [-1,1]
for i=1:n_y
    % X_norm 138*155
     
    input_scale(:,i)=((input(:,i)-min(input(:,i)))/(max(input(:,i))-min(input(:,i))))*2-1;
end

input_scale=input_scale;

% for i=1:n_y
% [p,h] = ranksum(input(:,i),label(:,1));
% P(i,:)=p;
% H(i,:)=h;
% end


for i=1:n_x
    if label(i)==1
        input_b(i_b,:)=input_scale(i,:);
        i_b=i_b+1;
    else %label(i)
        input_m(i_m,:)=input_scale(i,:);
        i_m=i_m+1;
    end
end

for i=1:n_y
[p,h] = ranksum(input_b(:,i),input_m(:,i));
% [h,p] = ttest2(input_b(:,i),input_m(:,i));
% [h,p] = vartest2(input_b(:,i),input_m(:,i));
P(i,:)=p;
H(i,:)=h;
end

clear Result;
for k=1 %p值的小数点位数
% hh = zeros(16384,1);
% hh(find(P<0.05)) = 1;
% H = hh;
% TempP=zeros(n_y,1);
j=1;
input_scale=input_scale';
for i=1:n_y
%     if(H(i)==1)
if (P(i,:)<0.01)
clear ouput_test
% if (P(i,:)<0.1^k)
        output_test(j,:)=input_scale(i,:);
%         TempP(i,:)=1;
        j=j+1;
        H2(i,:)=1;
else
%         Feature_Ttest2_BC(i,:)=[];
        H2(i,:)=0;
    end
   
end
input_scale=input_scale';
output_test=output_test';

% feature = output_test;
% % save('P_feature_1230_4','feature')
[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,output_test,label);
Result(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
save('TI_FeatureSelect_P_ttest2_1-10_10flod_0.00001.mat','output_test','label','H2','P','Result')

end
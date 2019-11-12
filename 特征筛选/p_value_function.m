function [output_test,P,H2,p_index]=p_value_function(input,label,p_range,method)
% clear output_test
% clear input_b
% clear input_m
% clear input
% clear P
% clear H
% clear H2
output_test=[];
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
    switch method
        case 'ranksum' 
            [p,h] = ranksum(input_b(:,i),input_m(:,i));
        case 'ttest2'
            [h,p] = ttest2(input_b(:,i),input_m(:,i));
        case 'vartest2'
            [h,p] = vartest2(input_b(:,i),input_m(:,i));
    end
P(i,:)=p;
H(i,:)=h;
end

for k=1 %p值的小数点位数
% hh = zeros(16384,1);
% hh(find(P<0.05)) = 1;
% H = hh;
% TempP=zeros(n_y,1);
j=1;
input_scale=input_scale';
for i=1:n_y
%     if(H(i)==1)
if (P(i,:)<p_range)
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
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,output_test,label);
% Result(k,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% save('TI_FeatureSelect_P_ttest2_1-10_10flod_0.00001.mat','output_test','label','H2','P','Result')
p_index=(find(H2==1));
end
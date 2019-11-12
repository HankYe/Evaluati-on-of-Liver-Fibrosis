function [temp_Feature_selection,flip_index,Data_sort]=feature_selection(feature,label,method)
features=feature;
label=label-1; %0,1 
% method='SRCnew';% SRC SRCnew LASSO mrmr Fisher_score reliefF F_score GA_mRMR

[imgRow,imgCol]=size(features);

switch method
    case 'SRC'
        %src
        [src_A,flip_index,Data_sort] = src_calculate(features,label);
    case 'SRCnew'
        % SRCnew
        [flip_index,Data_sort]=SRCnew(features,label+1);
    case 'LASSO'
        %LASSO
        [LASSO_wieghts,flip_index,Data_sort]=LASSO_calculate(features,label);
    case 'mrmr'
        % % mrmr
        [flip_index,Data_sort] = mrmr(features,label); %flip_index,Data_Fscore_sort相对应
    case 'Fisher_score'
        % Fisher_score
        [Fisher_Weight,flip_index,Data_sort]= fsFisher(features,label) ;
    case 'reliefF'
        %reliefF 要求label是1，2
        [relieF_Weight,flip_index,Data_sort]= reliefF(features,label+1,100,10,0,imgCol) ;
    case 'F_score'
        % F_score
        [F_score,flip_index,Data_sort] = F_score_calculate(features,label+1); %flip_index,Data_Fscore_sort相对应 label要求1,2
    case 'GA_mRMR'
        [gamrmr_output,time,value,ind,n_fea,g_sel,C_sel,temp,fea_sel] = libsvm_GA_mRMR_feature_select(features, label, 50, 30, 44, 0.9, 0.1, 1, 0, 10);
%         save('Data_GAmrmr','gamrmr_output','time','value','ind','n_fea','g_sel','C_sel','temp','fea_sel','label','feature');
        flip_index=find(temp==1);
        Data_sort=features(:,flip_index);
end

temp_Feature_selection(1,1:imgCol)=0;
temp_Feature_selection(1,flip_index)= 1;
        
        
end

% 
% %十分法
% % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,features,label);
% % 
% % % % % 原始准确度
% % % % % % % % % % % tic
% % % % % % % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,features,label);
% % % % % % % Result_488(1,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]% save('Result_origin197','Result_origin197')
% % % % % % % % % % % 
% % % % % toc
% 
% %src
% [src_A,flip_index,Data_src_A_sort] = src_calculate(features,label);
% for i=20%0:15%imgCol
% Data=Data_src_A_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_SRC_best2(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_SRC_best','Result_SRC_best');
% end
% 
% % SRCnew
% [SRCnew_index,Data_SRCnew_sort]=SRCnew(features,label-1);
% for i=20:5:90%:100%1:imgCol
% Data=Data_SRCnew_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_SRCnew_TENTH(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_mrmr_best','Result_mrmr_best');
% end
%  
% %LASSO
% [LASSO_wieghts,flip_index,Data_LASSO_sort]=LASSO_calculate(features,label);
% for i=20:5:90%:100%1:imgCol
% Data=Data_LASSO_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_LASSO_TENTH(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_mrmr_best','Result_mrmr_best');
% end
% 
% 
% % % mrmr
% [mrmr_index,Data_mrmr_sort] = mrmr(features,label); %flip_index,Data_Fscore_sort相对应
% for i=50%:100%1:imgCol
% Data=Data_mrmr_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% % [auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_mrmr_best(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_mrmr_best','Result_mrmr_best');
% end
% % % 
% % Fisher_score
% [Fisher_Weight,flip_index,Data_Fisher_sort]= fsFisher(features,label) ;
% for i=50%1:imgCol
% Data=Data_Fisher_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% %[auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_Fisher_best(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_Fisher_best','Result_Fisher_best');
% end
% 
% %reliefF 要求label是1，2
% % features=gamrmr_output;
% [imgRow,imgCol]=size(features);
% [relieF_Weight,flip_index,Data_relieF_sort]= reliefF(features,label+1,100,10,0,imgCol) ;
% for i=50%1:imgCol
% Data=Data_relieF_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% %[auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_relieF_best(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_relieF_best','Result_relieF_best');
% end
%  
% % F_score
% [F_score,flip_index,Data_Fscore_sort] = F_score_calculate(features,label+1); %flip_index,Data_Fscore_sort相对应 label要求1,2
% for i=50%1:imgCol
% Data=Data_Fscore_sort(:,1:i);
% [auc,accuracy,sensitivity,specificity,ppv,npv,mcc]=Adaboost_or_SVM_TENTH(5,1,Data,label);
% %[auc,accuracy,sensitivity,specificity,ppv,npv,mcc,TEMP]=Adaboost_or_SVM_LOOCV(1,1,Data,label);
% Result_Fscore_best(i,:)=[auc,accuracy,sensitivity,specificity,ppv,npv,mcc]
% % save('Result_Fscore_best','Result_Fscore_best');
% end

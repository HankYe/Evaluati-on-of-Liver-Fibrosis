function [auc,acc,sens,spec,ppv,npv,mcc]=assess(Label_predict,Label_golden,Predict)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Label_predict:分类器预测标签，1,2,1*N
% % Label_golden：原始真正金标准值，1,2,1*N
% % Predict：分类器决策值,scores
% %
% % auc: area under curve ROC曲线下的面积
% % acc: accuracy 准确率
% % sens: sensitivity 灵敏度
% % spec: specificity 特异度
% % ppv: positive predictive value 阳性预测率
% % npv: negative predictive value阴性预测率
% % mcc: Matthew's correlation coefficient 相关系数 
% %
% % TP:表示分类正确：True Positive：本来是正样例，分类成正样例。
% % FP:False Positive ：本来是负样例，分类成正样例，通常叫误报。
% % FN:False Negative：本来是正样例，分类成负样例，通常叫漏报。
% % TN: True Negative：本来是负样例，分类成负样例。表示分类错误：
% %
% % by：刘桐桐 2016.05.05 liutongtongsh@hotmail.com
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test heart scale demo
% load('heart_scale_labelandscores.mat')
% [ auc,acc,sens,spec,ppv,npv,mcc]= assess(predict_label',heart_scale_label',dec_values);
%
% test heart scale 在assess中添加程序段  
%     Label_predict(find(Label_predict==1))=2;    
%     Label_predict(find(Label_predict==-1))=1;
%     Label_golden(find(Label_golden==1))=2;
%     Label_golden(find(Label_golden==-1))=1;
%
% heart scale result
% acc =    0.8667
% sens =    0.8083
% spec =    0.9133
% ppv =    0.8818
% npv =    0.8562
% mcc =    0.7298
% auc =    0.9304  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Label_predict(find(Label_predict==1))=0;
    Label_predict(find(Label_predict==2))=1;
    Label_golden(find(Label_golden==1))=0;
    Label_golden(find(Label_golden==2))=1;


    TP = sum(Label_predict&Label_golden); % real - m, test - m
    FP = sum(Label_predict&~Label_golden);% real - b, test - m
    FN = sum(~Label_predict&Label_golden);% real - m, test - b
    TN = sum(~Label_predict&~Label_golden); % real - b, test - b  
    %要求label为0,1
    e=1e-10;

    acc = (TP+TN)/(TP+TN+FP+FN+e) % accuracy 准确率
    sens = TP/(TP+FN+e) % sensitivity 灵敏度
    spec = TN/(TN+FP+e) % specificity 特异度
    ppv = TP/(TP+FP+e) % positive predictive value 阳性预测率
    npv = TN/(TN+FN+e) % negative predictive value阴性预测率
    mcc = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)+e) % mcc: Matthew's correlation coefficient 相关系数 

    scale_label=Label_golden';
    dec_values=Predict';
    
%方法一 好用的
%    auc = plot_roc( dec_values, scale_label ) %要求进去的scale_label是0,1

%方法2 好用的   
   auc = plot_roc_and_det( dec_values, scale_label ) %要求进去的scale_label是0,1

end
function [auc,acc,sens,spec,ppv,npv,mcc]=assess(Label_predict,Label_golden,Predict)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Label_predict:������Ԥ���ǩ��1,2,1*N
% % Label_golden��ԭʼ�������׼ֵ��1,2,1*N
% % Predict������������ֵ,scores
% %
% % auc: area under curve ROC�����µ����
% % acc: accuracy ׼ȷ��
% % sens: sensitivity ������
% % spec: specificity �����
% % ppv: positive predictive value ����Ԥ����
% % npv: negative predictive value����Ԥ����
% % mcc: Matthew's correlation coefficient ���ϵ�� 
% %
% % TP:��ʾ������ȷ��True Positive�����������������������������
% % FP:False Positive �������Ǹ��������������������ͨ�����󱨡�
% % FN:False Negative��������������������ɸ�������ͨ����©����
% % TN: True Negative�������Ǹ�����������ɸ���������ʾ�������
% %
% % by����ͩͩ 2016.05.05 liutongtongsh@hotmail.com
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test heart scale demo
% load('heart_scale_labelandscores.mat')
% [ auc,acc,sens,spec,ppv,npv,mcc]= assess(predict_label',heart_scale_label',dec_values);
%
% test heart scale ��assess����ӳ����  
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
    %Ҫ��labelΪ0,1
    e=1e-10;

    acc = (TP+TN)/(TP+TN+FP+FN+e) % accuracy ׼ȷ��
    sens = TP/(TP+FN+e) % sensitivity ������
    spec = TN/(TN+FP+e) % specificity �����
    ppv = TP/(TP+FP+e) % positive predictive value ����Ԥ����
    npv = TN/(TN+FN+e) % negative predictive value����Ԥ����
    mcc = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)+e) % mcc: Matthew's correlation coefficient ���ϵ�� 

    scale_label=Label_golden';
    dec_values=Predict';
    
%����һ ���õ�
%    auc = plot_roc( dec_values, scale_label ) %Ҫ���ȥ��scale_label��0,1

%����2 ���õ�   
   auc = plot_roc_and_det( dec_values, scale_label ) %Ҫ���ȥ��scale_label��0,1

end
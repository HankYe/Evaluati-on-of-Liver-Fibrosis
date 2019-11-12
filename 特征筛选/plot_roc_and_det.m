function [auc]=plot_roc_and_det( dec_values, scale_label )

    [X,Y,T,AUC] =perfcurve(scale_label,dec_values,'1');
    %[X,Y,T,AUC] = perfcurve(labels,scores,posclass,'param1', val1,'param2',val2,...)
    %labels:Ŀ���ǩ scores:����ֵ posclass�������ǩ
    %'param'��val���Զ���X��Y�����ֵ��������Կ�����������Ĭ���Ƕ���X��ΪFPR��Y��ΪTPR
    auc=AUC;
   hold on 
    plot([0,1],[0,1],'b--') % ��һ���Խ��ߣ���(0,0)��(1,1)
    plot(X,Y)
%     Z=1-X+Y-1;M2=find(Z==max(Z));text(X(M2),Y(M2),'*','color','r');
%     hold on
%     
%     dn = 8;   % ��ǵ�����ÿ20������һ��
%     for i=1:length(X)
%         if(mod(i, dn)==0)
%             plot(X(i), Y(i), '^');
%         end
%     end
    %xlabel('sensitivity'),ylabel('1-specificity'),title('ROC curve');%���ROC����
    
%     [X,Y] =perfcurve(scale_label,dec_values,'1','xCrit','FPR','yCrit','FNR');
%     %���¶������ֵ��'xCrit','FPR'��ʾ����X���ΪFPR=FP/(TN+FP),������ʵ����Ĭ��ֵ
%     %'yCrit','FNR'��ʾ����y���ΪFNR=FN/(TP+FN)
% 
%     figure,plot(X,Y),xlabel('fall'),ylabel('miss'),title('DET');%���DET����
end
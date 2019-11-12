function [auc]=plot_roc_and_det( dec_values, scale_label )

    [X,Y,T,AUC] =perfcurve(scale_label,dec_values,'1');
    %[X,Y,T,AUC] = perfcurve(labels,scores,posclass,'param1', val1,'param2',val2,...)
    %labels:目标标签 scores:决策值 posclass：正类标签
    %'param'和val可以定义X和Y的输出值，具体可以看函数帮助，默认是定义X轴为FPR，Y轴为TPR
    auc=AUC;
   hold on 
    plot([0,1],[0,1],'b--') % 画一条对角线，从(0,0)到(1,1)
    plot(X,Y)
%     Z=1-X+Y-1;M2=find(Z==max(Z));text(X(M2),Y(M2),'*','color','r');
%     hold on
%     
%     dn = 8;   % 标记点间隔，每20个点标记一次
%     for i=1:length(X)
%         if(mod(i, dn)==0)
%             plot(X(i), Y(i), '^');
%         end
%     end
    %xlabel('sensitivity'),ylabel('1-specificity'),title('ROC curve');%输出ROC曲线
    
%     [X,Y] =perfcurve(scale_label,dec_values,'1','xCrit','FPR','yCrit','FNR');
%     %重新定义输出值，'xCrit','FPR'表示定义X输出为FPR=FP/(TN+FP),这里其实就是默认值
%     %'yCrit','FNR'表示定义y输出为FNR=FN/(TP+FN)
% 
%     figure,plot(X,Y),xlabel('fall'),ylabel('miss'),title('DET');%输出DET曲线
end
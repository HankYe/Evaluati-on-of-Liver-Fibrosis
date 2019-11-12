% clc;
% clear all;
% close all;
% 
% load Target138.mat;
% load Breast138_morph11_auto_new.mat;
% load PC_LBPV138_8ori_6scale_r2p16.mat

% load ./BUS 138/GAresult_gen30_fun1_4; %GA_mRMR结果
% load ./BUS 138/GAresult_gen30_fun0_1; %GA结果
function [AUC,ACC,SENS,SPEC,PPV,NPV,MCC] = libsvm_LOOCV(feature_input, label_target,g_sel,C_sel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 x = [feature_input'];
% x = [Breast138_morph11_auto_new' ; PC_LBPV138_8ori_6scale_r2p16'];
% x = [PC_LBPV138_8ori_6sale_r2p16'];
% x = [Breast138_morph11_auto_new'];

[no dim] = size(x');  % database

% feature normalization [-1,1]
for i=1:dim
    x(i,:)=((x(i,:)-min(x(i,:)))/(max(x(i,:))-min(x(i,:))))*2-1;
end

%%%%%%%%%%%%%%%%%%%%
X = x';
T_scale = label_target*2-1;
%T_scale = T138*2-1;

train_data_labels = T_scale';
% fea_sel = find(temp==1); % 找到所选择特征的索引值
fea_sel = [1:dim];
X_scale = X(:,fea_sel);
% 
% %%%%%%%%%%%%%%%%%%%本来这段是运行的
% % % default setting
% C_sel = 1;
% g_sel = 1/dim; % 1/#features

% 
% [no,dim]=size(label_target);
% for i=1:no
%     if (label_target(i)==1)
%         label_target(i)=1;
%     else label_target(i)=0;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sens_predict, spes_predict, accuracy,sensitivity,specificity,ppv,npv,mcc ] = myclassifysvmROC_LOOCV_new_pro( X_scale,X_scale,label_target,C_sel,g_sel);
% [sens_predict, spes_predict, accuracy,sensitivity,specificity,ppv,npv,mcc ] = myclassifysvmROC_LOOCV_new_pro( X_scale,X_scale,T138,C_sel,g_sel);
ACC = accuracy;
SENS = sensitivity;
SPEC = specificity;
PPV = ppv;
NPV = npv;
MCC = mcc;
% AUC = 1+sum(0.5*(sens_predict(2:end)+sens_predict(1:end-1)).*(spes_predict(2:end) - spes_predict(1:end-1)));
v1 = sens_predict; v2 = 1-spes_predict;
AUC = -sum(0.5*(v1(2:end)+v1(1:end-1)).*(v2(2:end) - v2(1:end-1)));
sensROC = sens_predict;
spesROC = spes_predict;
sens_mRMR=mean(sensROC,1);
spes_mRMR=mean(spesROC,1);
% plot(1-spes_mRMR,sens_mRMR,'LineWidth',1,'Color','b');hold on;
% xlabel('1 - 特异度','fontsize',10);ylabel('敏感度','fontsize',10);

% h = legend('全部特征','GA','GA_mRMR',4,'fontsize',10);
% set(h,'Interpreter','none')
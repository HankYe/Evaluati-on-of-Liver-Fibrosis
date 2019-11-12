function [Fitness_value Obj_value Rank_value] = Fitness_function_pro_libsvm(  Nind, Chrom, X_norm, dim, T , K_fold, g_value, C_value, fea_rank)

% ���������ʼ��
Fitness_value = zeros(Nind,1); % ��Ӧ�Ⱥ���
Obj_value = zeros(Nind,1);  % Ŀ�꺯����������׼ȷ��
Rank_value = zeros(Nind,1); % mRMR����

Chrom_fea = Chrom(:,1:dim);
% sampletotal=size(X_norm,1);
% 
% M = find(T==1); % malignant index 
% B = find(T==0); % benign index 
% 
% n_malignant = length(M); % malignant
% n_benign = sampletotal- n_malignant; % benign
% 
% % Rearrange X, [69 malignant, 69 benign]
% Mindex= M(randperm(n_malignant));
% X_malignant = X_norm(Mindex,:);
% T_malignant = ones(1,n_malignant);
% 
% Bindex= B(randperm(n_benign));
% X_benign = X_norm(Bindex,:);
% T_benign = -ones(1,n_benign);
% 
% indices_m = crossvalind('Kfold',n_malignant,K_fold);
% indices_b = crossvalind('Kfold',n_benign,K_fold);

for n = 1:Nind
    Chrom_fea_ind = find(Chrom_fea(n,:) == 1);
    Chrom_fea_rank = fea_rank(Chrom_fea_ind);
    
    X_Chrom = X_norm(:,Chrom_fea_ind);
    T_Chrom = T';
    T_scale = T_Chrom*2-1;
    
    
%     SVMModel = fitcsvm(X_Chrom,T_scale,'Standardize',true,'KernelFunction','RBF',...
%         'KernelScale','auto');
%     CVSVMModel = crossval(SVMModel);
%     classLoss = kfoldLoss(CVSVMModel);
%     Acc=(1-classLoss)*100;
    
    % �ֱ����ÿ��Chrom�ķ���׼ȷ�ԣ��õ��������Obj����
%     cmd = [ '-t 2 ','-v ',num2str(K_fold), ' -c ',num2str(C_value(n)), ' -g ', num2str(g_value(n)) , ' -q'];
    cmd = [ '-t 2 ','-v ',num2str(K_fold), ' -c ',num2str(C_value(n)), ' -g ', num2str(g_value(n)) ];

    Acc = libsvmtrain(T_scale, X_Chrom, cmd);  
    
    Obj_value(n) = Acc;
    
    Rank_value(n) = sum(Chrom_fea_rank)/sum(1:dim); %�����ֵ��mRMR��rank�͹�һ��
%     Fitness_value(n) = 1+Obj_value(n)-Rank_value(n);
%     Fitness_value(n) = 0.8*Obj_value(n)+(1-0.8)/(1+Rank_value(n));
%     Fitness_value(n) = Obj_value(n)+1/(1+Rank_value(n));
    Fitness_value(n) = Obj_value(n)+1-Rank_value(n)/(1+Obj_value(n));
end


end


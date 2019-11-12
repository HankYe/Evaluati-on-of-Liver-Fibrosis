function Fitness_value = Fitness_function_libsvm( Nind, Chrom, X_norm, dim, T , K_fold, g_value, C_value)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Fitness_value = zeros(Nind,1);
Chrom_fea = Chrom(:,1:dim);
% sampletotal=size(X_norm,1);

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
    X_Chrom = X_norm(:,Chrom_fea_ind);
    T_Chrom = T';
    T_scale = T_Chrom*2-1;
    % 分别计算每个Chrom的分类准确性，得到结果存于Fitness变量
    
      
%     SVMModel = fitcsvm(X_Chrom,T_scale,'Standardize',true,'KernelFunction','RBF',...
%         'KernelScale','auto');
%     CVSVMModel = crossval(SVMModel);
%     classLoss = kfoldLoss(CVSVMModel)
%     Acc=(1-classLoss)*100;
    
%     cmd = [ '-t 2 ','-v ',num2str(K_fold), ' -c ',num2str(C_value(n)), ' -g ', num2str(g_value(n)),'-q'];
     cmd = [ '-t 2 ','-v ',num2str(K_fold), ' -c ',num2str(C_value(n)), ' -g ', num2str(g_value(n)) ];

    Acc = libsvmtrain(T_scale, X_Chrom, cmd); %20170726 change by ltt
% % %     不需要libsvm了而是使用系统自带的函数fitcsvm
    
    Fitness_value(n) = Acc;

end


end


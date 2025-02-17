function multiSVMStruct = multiSVMTrain(TrainData, nSampPerClass, nClass, C, gamma)
%function multiSVMStruct = multiSVMTrain(TrainData, nSampPerClass, nClass, C, gamma)
% 采用1对1投票策略将 SVM 推广至多类问题的训练过程，将多类SVM训练结果保存至multiSVMStruct中
%
% 输入:--TrainData:每行是一个样本人脸
%     --nClass:人数，即类别数
%     --nSampPerClass:nClass*1维的向量，记录每类的样本数目，如 nSampPerClass(iClass)
%     给出了第iClass类的样本数目
%     --C:错误代价系数，默认为 Inf
%     --gamma:径向基核函数的参数 gamma，默认值为1
%
% 输出:--multiSVMStruct:一个包含多类SVM训练结果的结构体

% 默认参数
if nargin < 4
    C = Inf;
    gamma = 1;
elseif nargin < 5
    gamma = 1;
end



%开始训练，需要计算每两类间的分类超平面，共(nClass-1)*nClass/2个
for ii=1:(nClass-1)
    for jj=(ii+1):nClass
        clear X;
        clear Y;
        startPosII = sum( nSampPerClass(1:ii-1) ) + 1;
        endPosII = startPosII + nSampPerClass(ii) - 1;
        X(1:nSampPerClass(ii), :) = TrainData(startPosII:endPosII, :);
            
        startPosJJ = sum( nSampPerClass(1:jj-1) ) + 1;
        endPosJJ = startPosJJ + nSampPerClass(jj) - 1;
        X(nSampPerClass(ii)+1:nSampPerClass(ii)+nSampPerClass(jj), :) = TrainData(startPosJJ:endPosJJ, :);
        
        
        % 设定两两分类时的类标签
        Y = ones(nSampPerClass(ii) + nSampPerClass(jj), 1);
        Y(nSampPerClass(ii)+1:nSampPerClass(ii)+nSampPerClass(jj)) = 0;
        
        % 第ii个人和第jj个人两两分类时的分类器结构信息
        CASVMStruct{ii}{jj}= svm_train( X, Y, 'Kernel_Function', @(X,Y) kfun_rbf(X,Y,gamma), 'boxconstraint', C );
     end
end

% 已学得的分类结果
multiSVMStruct.nClass = nClass;
multiSVMStruct.CASVMStruct = CASVMStruct;

% 保存参数
save('Mat/params.mat', 'C', 'gamma');




function multiAdaboostStruct = multiAdaboostTrain(TrainData, nSampPerClass,nClass,iteration)
%function multiAdaboostStruct = multiAdaboostTrain(TrainData, nSampPerClass)
% 采用1对1投票策略将 Adaboost 推广至多类问题的训练过程，将多类Adaboost训练结果保存至multiAdaboostStruct中
%
% 输入:--TrainData:每行是一个样本人脸
%     --nClass:人数，即类别数
%     --nSampPerClass:nClass*1维的向量，记录每类的样本数目，如 nSampPerClass(iClass)
%     给出了第iClass类的样本数目
%
% 输出:--multiAdaboostStruct:一个包含多类Adaboost训练结果的结构体
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
        Y = ones(nSampPerClass(ii) + nSampPerClass(jj), 1)+1;
        Y(nSampPerClass(ii)+1:nSampPerClass(ii)+nSampPerClass(jj)) = 1;
        
        % 第ii个人和第jj个人两两分类时的分类器结构信息
        AdaboostStruct{ii}{jj}=ADABOOST_tr(@threshold_tr,@threshold_te, X, Y, iteration);%进行iteration次迭代，iteration个弱分类器
     end
end

% 已学得的分类结果
multiAdaboostStruct.AdaboostStruct = AdaboostStruct;
multiAdaboostStruct.nClass = nClass;
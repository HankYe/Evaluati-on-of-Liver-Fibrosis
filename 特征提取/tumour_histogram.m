function feature = tumour_histogram(tumour)
tumour_h = reshape(tumour, numel(tumour), 1);%%将肿瘤图像拉伸为一维向量
k = 1;
for i = 1:length(tumour_h)
    if tumour_h(i) ~=0
        tumour_adjust(k) = tumour_h(i);
        k = k+1;
    end
end
%% 常规参数计算
tmax = max(tumour_adjust); %%最大值
tmin = min(tumour_adjust); %%最小值
tmean = mean(tumour_adjust); %%平均值
tmedia = median(tumour_adjust); %%中值
tenergy = sum(tumour_adjust.^2);
N = length(tumour_adjust);
tkurtosis = (sum((tumour_adjust-tmean).^4)/N)/(sqrt(sum((tumour_adjust-tmean).^2)/N)).^2;
tMAD = sum(abs(tumour_adjust-tmean))/N;
tRange = tmax - tmin;
tRMS = sqrt(tenergy/N);
tskewness = (sum((tumour_adjust-tmean).^3)/N)/(sqrt(sum((tumour_adjust-tmean).^2)/N)).^3;
tdeviation = (sum((tumour_adjust-tmean).^2)/(N-1))^(1/2);
tvariance = sum((tumour_adjust-tmean).^2)/(N-1);
%tb = boxplot(tumour_adjust); %% 显示箱状图
figure;
th = histogram(tumour_adjust,100); %% 显示直方图
hold on;
%% 进行高斯拟合
x = (th.BinLimits(1)+th.BinWidth/2:th.BinWidth:th.BinLimits(2))';
y = (th.Values)';
% fitobject = fit(x,y,'gauss1');
% plot(fitobject);
close;
%% 计算6个直方图特征
h_entropy = 0; %%熵
y = y./length(tumour_adjust);%%对直方图进行归一化
h_mean = sum(x.*y); %% 直方图均值
h_variance = sum((x-h_mean).^2.*y);%%方差
h_skewness = sum((x-h_mean).^3.*y)/(h_variance^(3/2));%%歪斜度
h_kurtosis = sum((x-h_mean).^4.*y)/(h_variance^2)-3;%%峰态
h_uniformity = sum(y.^2);%%均匀度
for j = 1:size(x,1)
    if y(j)~=0
        h_entropy = h_entropy + y(j)*log2(y(j));
    end
end
h_entropy = - h_entropy;
%% 输出
% feature= [tenergy;h_entropy;tkurtosis;tmax;tmean;tMAD;tmedia;tmin;tRange;tRMS;tskewness;tdeviation;h_uniformity;tvariance;fitobject.a1;fitobject.b1;fitobject.c1;h_mean;h_variance;h_skewness;h_kurtosis];
feature= [tenergy;h_entropy;tkurtosis;tmean;tMAD;tmedia;tRange;tRMS;tskewness;tdeviation;h_uniformity;tvariance;h_mean;h_variance;h_skewness;h_kurtosis];
end
%% 希望找到直方图的两个峰值，这部分没有深入讨论
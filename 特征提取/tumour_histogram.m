function feature = tumour_histogram(tumour)
tumour_h = reshape(tumour, numel(tumour), 1);%%������ͼ������Ϊһά����
k = 1;
for i = 1:length(tumour_h)
    if tumour_h(i) ~=0
        tumour_adjust(k) = tumour_h(i);
        k = k+1;
    end
end
%% �����������
tmax = max(tumour_adjust); %%���ֵ
tmin = min(tumour_adjust); %%��Сֵ
tmean = mean(tumour_adjust); %%ƽ��ֵ
tmedia = median(tumour_adjust); %%��ֵ
tenergy = sum(tumour_adjust.^2);
N = length(tumour_adjust);
tkurtosis = (sum((tumour_adjust-tmean).^4)/N)/(sqrt(sum((tumour_adjust-tmean).^2)/N)).^2;
tMAD = sum(abs(tumour_adjust-tmean))/N;
tRange = tmax - tmin;
tRMS = sqrt(tenergy/N);
tskewness = (sum((tumour_adjust-tmean).^3)/N)/(sqrt(sum((tumour_adjust-tmean).^2)/N)).^3;
tdeviation = (sum((tumour_adjust-tmean).^2)/(N-1))^(1/2);
tvariance = sum((tumour_adjust-tmean).^2)/(N-1);
%tb = boxplot(tumour_adjust); %% ��ʾ��״ͼ
figure;
th = histogram(tumour_adjust,100); %% ��ʾֱ��ͼ
hold on;
%% ���и�˹���
x = (th.BinLimits(1)+th.BinWidth/2:th.BinWidth:th.BinLimits(2))';
y = (th.Values)';
% fitobject = fit(x,y,'gauss1');
% plot(fitobject);
close;
%% ����6��ֱ��ͼ����
h_entropy = 0; %%��
y = y./length(tumour_adjust);%%��ֱ��ͼ���й�һ��
h_mean = sum(x.*y); %% ֱ��ͼ��ֵ
h_variance = sum((x-h_mean).^2.*y);%%����
h_skewness = sum((x-h_mean).^3.*y)/(h_variance^(3/2));%%��б��
h_kurtosis = sum((x-h_mean).^4.*y)/(h_variance^2)-3;%%��̬
h_uniformity = sum(y.^2);%%���ȶ�
for j = 1:size(x,1)
    if y(j)~=0
        h_entropy = h_entropy + y(j)*log2(y(j));
    end
end
h_entropy = - h_entropy;
%% ���
% feature= [tenergy;h_entropy;tkurtosis;tmax;tmean;tMAD;tmedia;tmin;tRange;tRMS;tskewness;tdeviation;h_uniformity;tvariance;fitobject.a1;fitobject.b1;fitobject.c1;h_mean;h_variance;h_skewness;h_kurtosis];
feature= [tenergy;h_entropy;tkurtosis;tmean;tMAD;tmedia;tRange;tRMS;tskewness;tdeviation;h_uniformity;tvariance;h_mean;h_variance;h_skewness;h_kurtosis];
end
%% ϣ���ҵ�ֱ��ͼ��������ֵ���ⲿ��û����������
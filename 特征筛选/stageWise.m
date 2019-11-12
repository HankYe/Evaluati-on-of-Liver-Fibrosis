function [ wResult ] = stageWise( x, y, eps, runtime)  
    [m,n] = size(x);%���ݼ��Ĵ�С  
    wResult = zeros(runtime, n);%���յĽ��  
    w = zeros(n,1);  
    wMax = zeros(n,1);  
    for i = 1:runtime  
        ws = w'%���ÿһ�μ��������Ȩ��  
        lowestError = inf;%������Сֵ  
        for j = 1:n  
            for sign = -1:2:1  
                wTest = w;%��ʼ��  
                wTest(j) = wTest(j)+eps*sign;%ֻ�ı�һά����  
                yTest = x*wTest;  
                %�����  
                rssE = rssError(y, yTest);  
                if rssE < lowestError%����ã����滻  
                    lowestError = rssE;  
                    wMax = wTest;  
                end  
            end  
        end  
        w = wMax;  
        wResult(i,:) = w;  
    end  
end  
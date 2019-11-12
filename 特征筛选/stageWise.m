function [ wResult ] = stageWise( x, y, eps, runtime)  
    [m,n] = size(x);%数据集的大小  
    wResult = zeros(runtime, n);%最终的结果  
    w = zeros(n,1);  
    wMax = zeros(n,1);  
    for i = 1:runtime  
        ws = w'%输出每一次计算出来的权重  
        lowestError = inf;%定义最小值  
        for j = 1:n  
            for sign = -1:2:1  
                wTest = w;%初始化  
                wTest(j) = wTest(j)+eps*sign;%只改变一维变量  
                yTest = x*wTest;  
                %求误差  
                rssE = rssError(y, yTest);  
                if rssE < lowestError%如果好，就替换  
                    lowestError = rssE;  
                    wMax = wTest;  
                end  
            end  
        end  
        w = wMax;  
        wResult(i,:) = w;  
    end  
end  
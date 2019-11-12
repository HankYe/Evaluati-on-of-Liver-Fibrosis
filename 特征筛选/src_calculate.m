function [x,flip_index,Data_src_sort] = src_calculate(Data,label);
method = 2;

dim = size(Data,2);
for i = 1:dim
    Data(:,i)=((Data(:,i)-min(Data(:,i)))/(max(Data(:,i))-min(Data(:,i))))*2-1;   
end


Lp            =  0


maxIter       =  300;    %外循环截止条件
stopE         =  0.0001; %外循环截止条件
stopP         =  1e-4; %外循环截止条件

winLength     =  1/5;%5,10,20


Tal           =  0.004;
VarNum        =  110;

[imgRow,imgCol]=size(Data);


% % [x,y,goodIndex,goodSN,Err]=selectVars(feature,label,Lp,maxIter,stopE,stopP,winLength,Tal,VarNum);
[x,y,goodIndex,goodSN,Err]=selectVars2(Data,label,Lp,maxIter,stopE,stopP,winLength,Tal,VarNum,method);
goodIndex = goodIndex';
output    = Data(:,goodIndex);
%  save('feature_sel_src2_300_5','goodIndex','output')


[src_sort,index]=sort(abs(x));
flip_index=flipud(index)'; %第i个为src第i高的特征的序号  由高到低

Data_src_sort=zeros(imgRow,imgCol);
for i=1:imgCol
    Data_src_sort(1:imgRow,i)= Data(1:imgRow,flip_index(i));
end

end

















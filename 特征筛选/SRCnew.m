% clear;
% clc;
function [SRCnew_index2,Data_SRCnew_sort]=SRCnew(input,label)
%要求：特征除去NAN 0， label 1,2，1在前，2在后
index_nan=find(~isnan(mean(input)));
index_zero=find(mean(input)~=0);

index_nanandzero=intersect(index_nan,index_zero) ;
features_nozeroandnan=input(:,index_nanandzero);
input=features_nozeroandnan;


[n_x,n_y]=size(input);

% Feature normalization [-1,1]
for i=1:n_y
    % X_norm 138*155
    input_scale(:,i)=((input(:,i)-min(input(:,i)))/(max(input(:,i))-min(input(:,i))))*2-1;
end
% load('label.mat');
gnd = label;
% clear label  % 标签[1 2....]
% load('feature_norm.mat');
fea = input_scale;    % 每个特征归一化到[-1 1]

index1=find(label==1);
index2=find(label==2);

gnd(1:length(index1))=1;
gnd(  (length(index1)+1):(length(index1)+length(index2) ) )=2;

for i=1:length(index1)
fea1(i,:)=input(index1(i),:);
end

for i=1:length(index2)
fea2(i,:)=input(index2(i),:);
end
fea=[fea1;fea2];

SRCnew_index = Untitled6(gnd,fea);

SRCnew_index2=index_nanandzero(SRCnew_index);

for i=1:n_y
Data_SRCnew_sort(:,i)=input(:,SRCnew_index(i));
end
end
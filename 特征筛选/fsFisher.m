function [W,index,Data_Fisher_sort] = fsFisher(Data,Label)
%Fisher Score, use the N var formulation
%   X, the data, each raw is an instance
%   Y, the label in 1 2 3 ... format
 
numClass = max(Label);
[numData, numFeature] = size(Data);
out.W = zeros(1,numFeature);

% statistic for classes
cIDX = cell(numClass,1);
n_i = zeros(numClass,1);
for j = 1:numClass
    cIDX{j} = find(Label(:)==j);
    n_i(j) = length(cIDX{j});
end

% calculate score for each features
for i = 1:numFeature
    temp1 = 0;
    temp2 = 0;
    f_i = Data(:,i);
    u_i = mean(f_i);
    
    for j = 1:numClass
        u_cj = mean(f_i(cIDX{j}));
        var_cj = var(f_i(cIDX{j}),1);
        temp1 = temp1 + n_i(j) * (u_cj-u_i)^2;
        temp2 = temp2 + n_i(j) * var_cj;
    end
    
    if temp1 == 0
        out.W(i) = 0;
    else
        if temp2 == 0
            out.W(i) = 100;
        else
            out.W(i) = temp1/temp2;
        end
    end
end

[~, out.fList] = sort(out.W, 'descend');
out.prf = 1;

W=out.W;
index=(out.fList)';

Data_relieF_sort=zeros(numData,numFeature);
for i=1:numFeature
    Data_Fisher_sort(1:numData,i)= Data(1:numData,index(i));
end


end

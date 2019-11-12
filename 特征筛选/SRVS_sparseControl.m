function [goodIndex,goodSN,Err]=SRVS_sparseControl(A,X,y,tal,upLimit,beDisp)
%this function selected the columns according to the contributionFactor
%A is the original data with the last row is the order of the data 
%X is the solution of Ax=y with min Lp(x)
%contributionFactor is defines the contribute of a column of A to the
%construction of y with aj*xj
%Tal=norm(y-Ax)/norm(y);
%upLimit is the uplimit of the number of variables to be selected
% beDisp spesifies if to display the Error term against the number of
% variables selected
%goodIndex is the index extracted from the last row of A

if nargin<5
    upLimit=inf;
end
if nargin<6
    beDisp=0;
end
if nargin<4
    Alpha=0;
end


absX=abs(X);
[absX,indexX]=sort(absX,'descend');
indexGood=find(absX>0);
lindexGood=length(indexGood);
lX=length(X);
x0=zeros(lX,1);
Err=zeros(lindexGood,1);

goodSN=0;
for i=1:lindexGood
    x0(indexX(1:i))=X(indexX(1:i));
    Err(i)=norm(y-A*x0)/norm(y);
    if Err(i)<tal
        goodSN=i;
        break;
    else
        goodSN=lindexGood;
    end
end

if upLimit<goodSN
    goodSN=upLimit;
end

goodIndex=indexX(1:goodSN);
Err=Err(1:goodSN);

if beDisp==1
    figure;
    plot(Err);
    xlabel('Number of variables selected');
    ylabel('Error term coefficients(Tal)');
end



    

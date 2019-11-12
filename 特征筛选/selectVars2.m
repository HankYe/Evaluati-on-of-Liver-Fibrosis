function [xk1,y,goodIndex,goodSN,Dk]=selectVars2(data,y,Lp,maxIter,stopE,stopP,winLength,Tal,VarNum,method)
%this function perfroms the variable selection using the SRVS algorithm
%data is in the size of m*n, where m is the number of samples, n is the
%number of variables;
%y is in the size of n*1
%Lp specify the penalization term: using Lp norm
%maxIter defines the maximum iteration numbers; There gone be
%maxIter/floor(winLength) number of sub-matrix being extracted and solved
%StopE is the step difference threshold; 
%stopP is threshold for the probabilty that any given pair of column
%variables in data will be compared >(1-stopP) 
%winLength=n/n1; where n1 is the length of the sub-matrix length
%Tal=norm(y-Ax)/norm(y) is used for sparsity control
%VarNum is the number of variables to be selected
%nargin is the number of variables inputed
if nargin<9,
    VarNum=inf;
end
if nargin<8,
    Tal=0.4;
end
if nargin<7,
    winLength=1/20;
end
if nargin<6,
    stopP=1e-4;
end
if nargin<5,
    stopE=1e-3;
end
if nargin<4
    maxIter=1/winLength*3;
end
if nargin<3,
    Lp=0.5;
end


% OStatus=OStatus
[SN,DL]=size(data);

windL=floor(SN*winLength);%the length of the window
nTimes=0;

nu(1:SN) = 1:SN;
for i = 1:maxIter
    ind = (FisherYatesShuffle(nu))';
    index(1+(i-1)*SN:i*SN,1) = ind;
end

x=zeros(DL,1);%initial x for Ax=y
xk1=x;

stepL=windL;%the step length for the window moving

stepN=1;%the number of moving steps
while 1  
    if maxIter<stepN
        break;
    end
    xk0 = xk1;
    winStart=stepL*(stepN-1)+1;
    winEnd=winStart+windL-1;    
    if winEnd>length(index)
        break;
    end
    
    id = index(winStart:winEnd);
    A = data(id,:);  
    Y = y(id);
%从输入data中随机选取1/5作为数据输入
    if Lp==1
        stepX =Homotophy(A,Y,DL);
    elseif Lp==0
        stepX= SolveOMP(A, Y,DL);
    else
        stepX=myMFOCUSS(A, Y, Lp);
    end
    coe(:,stepN) = stepX;
    xk1 = mean(coe,2);
    dk = norm(xk1-xk0);
    
    
    if  dk<stopE  
        break;
    end
    
    Dk(stepN) = dk; 
    stepN=stepN+1;
end

% z = 1:size(Dk,2);
% figure
% plot(z,Dk)
% xlabel('iteration k')
% ylabel('d(k)')
% axis([0 350 0 0.5])

if  method == 1
        [s d] = sort(abs(xk1),'descend');
        goodIndex = d(1:VarNum);
        goodSN    = VarNum;
else
        if nargout>2
            %analyze the results
            [goodIndex,goodSN,Err]=SRVS_sparseControl(data,xk1,y,Tal,VarNum,0);  %最后一位是0则不显示曲线
        end
end
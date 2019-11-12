function [x,y,goodIndex,goodSN,Err]=selectVars(data,y,Lp,maxIter,stopE,stopP,winLength,Tal,VarNum)
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

data(SN+1,:)=1:DL;%the index of the original data
windL=floor(DL*winLength);%the length of the window
nTimes=0;

x=zeros(DL,1);%initial x for Ax=y
% y=zeros(SN,1);%set y for Ax=y
% y(1:class1Num)=1;%SN=[class1Num+class2Num]
x1=x;
% stopE=0.001;
% stopP=1e-4;
% Lp=0.5;
nIter=0;
while 1
    
    nIter=nIter+1;
    if maxIter<nIter
        break;
    end
    
    close all;
    x0=x1;%record the original one    
    A=FisherYatesShuffle(data);
    [stepx,stepErr]=stepSRVS_selectVars(A,windL,y,Lp);
    x=x+stepx;
    disp('------------------------');
    nTimes=nTimes+1;
    disp(sprintf('Iteration %d',nTimes));
    disp(sprintf('Step error=%2.2f',stepErr));
    x1=x/nTimes;
    stepDif(nTimes)=norm(x1-x0);
%     pNotEvalued(nTimes)=((((DL+1)/windL)-2)/(((DL+1)/windL)-1))^floor(nTimes); % the probability that a column vector is not elvaluated    
% probability that any two vector pairs are not compared
    pCompared(nTimes)=1-(1-winLength)^nTimes;
    disp(sprintf('Step difference=%2.5f,P_compared=%2.5f',stepDif(nTimes),pCompared(nTimes)));
    
    if stepDif(nTimes)<stopE && pCompared(nTimes)>(1-stopP)  
        break;
    end
        
end
disp('------------------------');
disp('The iteration stops!');
x=x/nTimes;

%test the solution
A=data(1:SN,:);
TotalErr=norm(A*x-y)

if nargout>2
    %analyze the results
    [goodIndex,goodSN,Err]=SRVS_sparseControl(A,x,y,Tal,VarNum,1);
end
function [x,stepErr]=stepSRVS_selectVars(A0,windL,y,Lp)
%windL is the number of columns in each windows
%class1Num is the number of samples in the first class, should corespond to
%the first class1Num rows
%A is the data to be processed, the last row of A0 contains the
%index each column of A
%Lp specifies the norm of the penalty term 

[rA,cA]=size(A0);
stepL=windL;%the step length for the window moving
x=zeros(cA,1);;%initial x for Ax=y
% y=zeros(rA-1,1);%set y for Ax=y
% y(1:class1Num)=1;%rA=[class1Num+class2Num]

stepN=1;%the number of moving steps

while 1   
    winStart=stepL*(stepN-1)+1;
    winEnd=winStart+windL-1;    
    if winEnd>cA
        break;
    end
    A=A0(1:rA-1,winStart:winEnd);  
    if Lp==1
        stepX =Homotophy(A,y,stepL);
    elseif Lp==0
        stepX= SolveOMP(A, y,stepL);
    else
        stepX=myMFOCUSS(A, y, Lp);
    end
    x(winStart:winEnd)=stepX;
    stepN=stepN+1;
end
x=x/(stepN);
stepErr=norm(A0(1:rA-1,:)*x-y);


oldOrder=A0(rA,:);%the original order
[oldOrder,indexOld]=sort(oldOrder,'ascend');
x=x(indexOld);%the solution of x






%use sparse solution to get the 
% [sx,numIters] =Homotophy(A,y,lx);
% [sx]= SolveBP(A, y,lx);
% [sols, numIters] = SolveLasso(A,y);
% [sols, numIters]= SolveStOMP(A, y,FLN*nClass);
% [sx, numIters]= SolveOMP(A, y,lx);
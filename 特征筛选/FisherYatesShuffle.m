function [A1,indexA]=FisherYatesShuffle(A)
%this function randomize the order of the column vectors of A according to
%Fisher-Yates shuffle algorithm
%indexA is the index after sha

[Ar,Ac]=size(A);

indexA=1:Ac;
for j=1:Ac
    jr=getRand(Ac);
    temp=indexA(j);
    indexA(j)=indexA(jr);
    indexA(jr)=temp;   
end
A1=A(:,indexA);


function ra=getRand(b)
% This function get a randomized int value between[1,b]
ra=rand(1,1);
ra=floor(ra*(b-1)+1);
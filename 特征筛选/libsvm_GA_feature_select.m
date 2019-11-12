function [output,time,value,ind,n_fea,g_sel,C_sel,temp,fea_sel] = libsvm_GA_feature_select(feature_input, label_target, Nind, Maxgen, Preci_add, Rxov, Rmut, GGAP, gen, K_fold)
% test for genetic algorithm optimization
% example:
% Nind = 50; %种群规模/个体数目
% Maxgen = 30; %遗传代数
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Preci = dim+24; % 变量的二进制位数，91fea+12sigma+12C
% Preci = dim+44; % 变量的二进制位数，91fea+12sigma+12C
% Rxov = 0.9; % 交叉概率
% Rmut = 0.1;% 变异概率
% GGAP = 1; % generation gap
% gen = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 X = [feature_input']';
% % X = [Breast138_morph11_auto_new' ; PC_LBPV138_8ori_6scale']';
%  X = [Breast138_morph11_auto_new' ; PC_LBPV138_8ori_6scale_r2p16']';

%%%原始的
% % X = [Breast138_morph11_auto_new' ; PC_LBPV138_8ori_6scale']';
% X = [Breast138_morph11_auto_new' ; PC_LBPV138_8ori_6scale_r2p16']'; 
 
[no,dim] = size(X);

% Feature normalization [-1,1]
for i=1:dim
    % X_norm 138*155
    X_norm(:,i)=((X(:,i)-min(X(:,i)))/(max(X(:,i))-min(X(:,i))))*2-1;
end


% GA optimization
Nind = Nind; %种群规模/个体数目
Maxgen = Maxgen; %遗传代数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Preci = dim+24; % 变量的二进制位数，91fea+12sigma+12C
Preci = dim+Preci_add; % 变量的二进制位数，91fea+12sigma+12C
Rxov = Rxov; % 交叉概率
Rmut = Rmut;% 变异概率
GGAP = GGAP; % generation gap
gen = gen; 
% % GA optimization
% Nind = 50; %种群规模/个体数目
% Maxgen = 30; %遗传代数
% %%%Preci = dim+24; % 变量的二进制位数，91fea+12sigma+12C
% Preci = dim+44; % 变量的二进制位数，91fea+12sigma+12C
% Rxov = 0.9; % 交叉概率
% Rmut = 0.1;% 变异概率
% GGAP = 1; % generation gap
% gen = 0;
% 生成初始种群
BaseV = crtbase(Preci,2);
[Chrom,Lind,BaseV] = crtbp(Nind,BaseV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 计算SVM中gamma及C的实数值，在Chrom中的后Preci_add个bin分别用于表示gamma和C
FieldD_g = [Preci_add/2;0.001;1;0;0;1;1]; % gamma:0.001-1
FieldD_C = [Preci_add/2;1;3000;0;0;1;1]; % C: 1-30,000
% % 计算SVM中gamma及C的实数值，在Chrom中的后24个bin分别用于表示gamma和C
% FieldD_g = [12;0.001;1;0;0;1;1]; % gamma:0.001-1
% FieldD_C = [12;1;3000;0;0;1;1]; % C: 1-30,000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始种群中的sigma值
g_value = bs2rv(Chrom(:,Preci-(Preci_add-1):Preci-Preci_add/2),FieldD_g);
C_value = bs2rv(Chrom(:,Preci-(Preci_add/2-1):end),FieldD_C);
% g_value = bs2rv(Chrom(:,Preci-23:Preci-12),FieldD_g);
% C_value = bs2rv(Chrom(:,Preci-11:end),FieldD_C);



% 使用k-fold cross validation计算分类的accuracy作为GA的fitness function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = label_target;
%T = T138;
K_fold = K_fold;%%%=10

% 初始种群的Fitnesss function calculation
FitV = Fitness_function_libsvm(  Nind, Chrom, X_norm, dim, T , K_fold, g_value, C_value);
tic
while gen < Maxgen

    % Selection
    SelCh = select('sus', Chrom, FitV, GGAP);
    
    % Crossover / Recombination
    SelCh = recombin('xovdp', SelCh, Rxov);
    
    % Mutation
    SelCh = mut(SelCh,Rmut);
    
    % Evaluation of offspring
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g_value_SelCh = bs2rv(SelCh(:,Preci-(Preci_add-1):Preci-Preci_add/2),FieldD_g);
    C_value_SelCh = bs2rv(SelCh(:,Preci-(Preci_add/2-1):end),FieldD_C);
%     g_value_SelCh = bs2rv(SelCh(:,Preci-23:Preci-12),FieldD_g);
%     C_value_SelCh = bs2rv(SelCh(:,Preci-11:end),FieldD_C);
    FitVSel = Fitness_function_libsvm( Nind*GGAP, SelCh, X_norm, dim, T, K_fold, g_value_SelCh, C_value_SelCh );
    
    % Reinsert offspring into population
    [Chrom FitV] = reins(Chrom,SelCh,1,1,FitV,FitVSel);
    
    % 遗传代数加1
    gen = gen+1;
    
end
time = toc
[value ind] = max(FitV)
temp = Chrom(ind,1:dim);
n_fea = sum(temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_bin = Chrom(ind,dim+1:dim+Preci_add/2);
C_bin = Chrom(ind,dim+Preci_add/2+1:end);
% g_bin = Chrom(ind,dim+1:dim+12);
% C_bin = Chrom(ind,dim+13:end);
g_sel = bs2rv(g_bin,FieldD_g)
C_sel = bs2rv(C_bin,FieldD_C)
fea_sel = find(temp==1);%%%选出特征的序号位置
temp;
input=feature_input;
ii=1;
for i=1:dim
    if (temp(i)==1)
        output(1:no,ii)=input(1:no,i);
        ii=ii+1;    
    else
    end
end
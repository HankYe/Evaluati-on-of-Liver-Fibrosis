warning off;
% clear

% load('TI_ratio_N5_87.mat')

currentpath = cd;

matpath = [currentpath, '\分割_非PTC\'];   
matfiles = dir(matpath);
file_no = size(matfiles,1);

addpath(genpath('radiomics-master\TextureToolbox'));

for i =3:file_no

    filename = matfiles(i).name;
    [path,name] = fileparts(filename);
    load([matpath,filename]); %

    BW3=BW3_tumor;BW3=logical(BW3);
    BW3_elastic=logical(BW3_elastic);
    BW3=logical(BW3);

    [FeaturesBase]=twoside_GetTextureFeatures(Img,BW3);
    
    %% 正式计算所有特征
    [ImgA1,ImgH1,ImgV1,ImgD1]=dwt2(Img,'coif1'); %进行小波变换,返回系数矩阵
    BW3_1=imresize(BW3,size(ImgA1));
 
    TextureHistFeaturesA_1=twoside_GetTextureFeatures(ImgA1,BW3_1);
    TextureHistFeaturesH_1=twoside_GetTextureFeatures(ImgH1,BW3_1);
    TextureHistFeaturesV_1=twoside_GetTextureFeatures(ImgV1,BW3_1);    
    TextureHistFeaturesD_1=twoside_GetTextureFeatures(ImgD1,BW3_1);

    FeatureAll=[FeaturesBase TextureHistFeaturesA_1 TextureHistFeaturesH_1 TextureHistFeaturesV_1 TextureHistFeaturesD_1   ];
%    save([currentpath,'\特征_PTCsecond\FeatureOf\',num2str(name,'%03d'),'.mat'],'FeatureAll','FeaturesBase','TextureHistFeaturesA_1','TextureHistFeaturesH_1','TextureHistFeaturesV_1','TextureHistFeaturesD_1','BW3','Img');
    Feature_TI_System_elastic(i-2,:)=[FeaturesBase TextureHistFeaturesA_1(:,1:end) TextureHistFeaturesH_1(:,1:end) TextureHistFeaturesV_1(:,1:end) TextureHistFeaturesD_1(:,1:end)];
end

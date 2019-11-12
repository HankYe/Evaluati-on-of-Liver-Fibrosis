function [TextureHistFeatures]=twoside_GetTextureFeatures(Img,BW3)
%     figplot=getfigplot(BW3);
    tumour=double(Img).*double(BW3);
    feature_histogram = double(tumour_histogram(tumour));
    figplot=getfigplot2(BW3);
    
 %%
    tumor=double(Img).*double(BW3);  %肿瘤区域提取  
    
    glcm = graycomatrix(tumor,'G',[]);
    out=GLCM_Features1(glcm,0); 
    glcmfeatures=[out.autoc,out.contr,out.corrm,out.corrp,out.cprom,out.cshad,out.dissi,out.energ,out.entro,out.homom,out.homop,...
       out.maxpr,out.sosvh,out.savgh,out.svarh,out.senth,out.dvarh,out.denth,out.inf1h,out.inf2h,out.homom,out.indnc,out.idmnc];
  %%
   [ROIonly,levels] = prepareVolume(Img,BW3,'other',1,1,1,1,'Matrix','Uniform',32);   
%     [GLCM] = getGLCM(ROIonly,levels); 
%     [glcmTextures] = getGLCMtextures(GLCM);
    [GLRLM] = getGLRLM(ROIonly,levels); 
    [glrlmTextures] = getGLRLMtextures(GLRLM); %输入为Gray-Level Run-Length Matrix灰度级运行矩阵，输出为纹理
    [GLSZM] = getGLSZM(ROIonly,levels); 
    [glszmTextures] = getGLSZMtextures(GLSZM);
    [NGTDM,countValid] = getNGTDM(ROIonly,levels); 
    [ngtdmTextures] = getNGTDMtextures(NGTDM,countValid);
    
   %%  
    TextureHistFeatures=double([glcmfeatures,(struct2array(glrlmTextures)),(struct2array(glszmTextures)),(struct2array(ngtdmTextures))...
         feature_histogram']);
   
end
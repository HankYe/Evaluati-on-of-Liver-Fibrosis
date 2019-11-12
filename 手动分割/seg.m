% clear all;
close all;



% clear
currentpath = cd;
matpath = [currentpath, '\123\']; % 自动
matfiles = dir(matpath);
file_no = size(matfiles,1);
for i =3 %(2+1):file_no% 3：081 +2？ 4：2+175
    flag=1;
    % for i=9 47 62 68 72 84   77
    
    filename = matfiles(i).name;
    [path,name] = fileparts(filename);
    %     Img=dicomread([matpath,filename]); %
    
% %     %读取分好的数据
%      path_read=[currentpath,'\1234\', 'seg', filename,'.mat'];
%     load(path_read);
    
    %%% 打开文件,ROI选择 %%%
    metadata = dicominfo([matpath,filename]);
%     y1=metadata.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0;y2=metadata.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1;x1=metadata.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0;x2=metadata.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1;
%     
    ImgDicom = dicomread([matpath,filename]);
    Img=ImgDicom;
    if size(Img,3)==1
    else
        Img=rgb2gray(Img);
    end
    ImgRaw=Img;
    path_save=[currentpath,'\1234\', 'seg' filename];
%     %             path_save=[currentpath,'\甲状腺随访第四批选出20170517分割结果\', filename];
    
%     if (flag==1)
%         %添一块
%         BW4= roipoly(Img);
%         BW5=BW3|BW4;BW3=BW5;
%         
%         nameS=[path_save  '.mat'];
%         %     save(nameS,'ImgDicom','Img','ImgRaw','y1','y2','x1','x2');
%         
%         save(nameS,'BW3','ImgDicom','Img','ImgRaw','y1','y2','x1','x2');
%         figure(1);imshow(Img,[]);
%         hold on;contour(BW3,[0 0],'.r','linewidth',1);
%         name=[ path_save   '.jpg'];
%         saveas(figure(1),name);
%         
%     else if(flag==2)
%             %去一块
%             BW4= roipoly(Img);
%             BW6=BW3&BW4;
%             BW5=BW3-BW6;
%             BW3=BW5;
%             
%             nameS=[path_save  '.mat'];
%             %     save(nameS,'ImgDicom','Img','ImgRaw','y1','y2','x1','x2');
%             
%             save(nameS,'BW3','ImgDicom','Img','ImgRaw','y1','y2','x1','x2');
%             figure(1);imshow(Img,[]);
%             hold on;contour(BW3,[0 0],'.r','linewidth',1);
%             name=[ path_save   '.jpg'];
%             saveas(figure(1),name);
%         else
            
            %保存原图
            % %         BW3= roipoly(Img);
            %
            % %     nameS=[path_save  '.mat'];
            % %     save(nameS,'BW3','ImgDicom','Img','ImgRaw');
            %     figure(1);imshow(Img,[]);
            %     hold on;contour(BW3,[0 0],'.r','linewidth',1);
            %     name=[ path_save   '.jpg'];
            %     saveas(figure(1),name);
            % % %     imshow(BW3)
            
            
            BW3= roipoly(Img);
            nameS=[path_save  '.mat'];
            %     save(nameS,'ImgDicom','Img','ImgRaw','y1','y2','x1','x2');
            
            save(nameS,'BW3','ImgDicom','Img','ImgRaw');
            figure(1);imshow(Img,[]);
            hold on;contour(BW3,[0 0],'.r','linewidth',1);
            %             name=[ path_save   '.jpg'];
            %             saveas(figure(1),name);

%             if (label==0)
%                 name=[ currentpath,'\1234\', 'seg' filename  '\不转移\'  '.jpg'];
%                 saveas(figure(1),name);
%             else
%                 name=[  path_save  '\转移\'  '.jpg'];
%                 saveas(figure(1),name);
%             end
            %     imshow(BW3)
            
        end
%     end
        close all
%     end
    
    
    %
    %
    % % imagepath='G:\Data_1st\US images both bmp&DICOM selected\data_select\1_2';
    % % Img= dicomread(imagepath );
    % imagepath='E:\科研\breast cancer\2.图像分割\手动分割\9.bmp';
    % Img= imread(imagepath );
    %
    % %% segment by hand
    % % path_save='G:\Data_1st\US images both bmp&DICOM selected\data_select\segment\img1_2';
    % path_save='E:\科研\breast cancer\2.图像分割\手动分割\img1_2';
    % BW= roipoly(Img);
    % nameS=[path_save  '.mat'];
    % save(nameS,'BW');
    % figure(1);imshow(Img,[]);hold on;contour(BW,[0 0],'.w','linewidth',1);
    % name=[path_save  '.tif'];
    % saveas(figure(1),name);
    %
    % %  load('G:\Data_1st\US images both bmp&DICOM selected\data_select\segment\img1_2.mat')
    % %  load('E:\科研\breast cancer\2.图像分割\手动分割\1.mat')
    % imshow(BW)
    
    

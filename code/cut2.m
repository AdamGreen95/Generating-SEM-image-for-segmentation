%%
close all
cd /home/hadmatter/桌面/homogeneous_nanoparitcle_cluster/imagedata/dataset
% folder strcture
% dataset
%     >|->train
%         |-data_0
%             |-img
%                 |-*.png
%             |-mask
%                 |-*.png
%                 |-*.png
%         |-img_1
%         ....
%     >|->test
%         |-*.png

%%
% folder dataset
clearvars -except map map_mask
dirOutput=dir('*.mat');
fileNames={dirOutput.name}';
filenum = numel(fileNames);
% save the varibles to mat file
filename = [num2str(filenum),'.mat']
save(filename,'map','map_mask');

%%
% create folder to store one image
cd train
%
basicsize = 256;
cutp1 = [0 0 basicsize basicsize];
cutp2 = [250 0 basicsize basicsize];
cutp3 = [472 0 basicsize basicsize];
cutp4 = [0 250 basicsize basicsize];
cutp5 = [250 250 basicsize basicsize];
cutp6 = [472 250 basicsize basicsize];
cutp = [cutp1;cutp2;cutp3;cutp4;cutp5;cutp6];

%%
step = length(cutp);
startpoint = filenum*step;
%
for i = 1:step
    trnfoldername = ['data_',num2str(startpoint)];
    mkdir(trnfoldername)
    cd(trnfoldername)
    
    %save the cropped image patch
    imagetemp = imcrop(map,cutp(i,:));
    mkdir img
    cd img
    imname = ['image_',num2str(filenum),'_',num2str(i),'.png'];
    imwrite(imagetemp,imname);
    cd ..
    
    %save its cooresponding masks into one folder
    masktemp = imcrop(map_mask,cutp(i,:));
    imzero = double(zeros(size(masktemp)));
    ptcnum = unique(masktemp(:));
    mkdir mask
    cd mask
    for j = 1:length(ptcnum)
        if j~=1
            imzero2 = imzero;
            imzero2(find(masktemp==ptcnum(j)))=1;
            mskname = ['image_',num2str(filenum),'_',num2str(i),...
                '_mask_',num2str(j),'.png'];
            imwrite(imzero2,mskname);
        end
    end
    cd ..
    startpoint = startpoint+1;
    cd ..
end


cd /home/hadmatter/桌面/homogeneous_nanoparitcle_cluster







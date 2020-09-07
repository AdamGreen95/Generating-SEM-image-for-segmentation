%%
close all
cd
cd imagedata/image
dirOutput=dir('*.png');
fileNames={dirOutput.name}'
filenum = numel(fileNames);
cd ..
%%

basicsize = 256;
cutp1 = [0 0 basicsize basicsize];
cutp2 = [250 0 basicsize basicsize];
cutp3 = [472 0 basicsize basicsize];
cutp4 = [0 250 basicsize basicsize];
cutp5 = [250 250 basicsize basicsize];
cutp6 = [472 250 basicsize basicsize];
cutp = [cutp1;cutp2;cutp3;cutp4;cutp5;cutp6];

%%

[maskx,masky] = gradient(double(map_mask));
mask = zeros(size(map));
mask(find(maskx~=0))=1;
mask(find(masky~=0))=1;
imshowpair(map,mask);

%%

for i = 1:6
    imagetemp = imcrop(map,cutp(i,:));
    masktemp = imcrop(mask,cutp(i,:));
    if filenum==0
        imagename = ['0.png'];
    else
        imagename = [num2str(filenum),'.png']
    end
    cd image
    imwrite(imagetemp,imagename);
    cd ..
    cd label
    imwrite(masktemp,imagename);
    cd ..
    filenum = filenum+1; 
end
cd ..


%%
cd /home/hadmatter/文档/PycharmProjects/unet-master/data/particle/test
dirim = dir('*.png');
num = numel(dirim)/2;
figure;
for i = 1:num
    imname = [num2str(i-1),'.png'];
    prename = [num2str(i-1),'_predict.png'];
    im = imread(imname);
    pre = imread(prename);
    subplot(2,3,i)
    imshowpair(im,pre)
    %imshowpair(im,pre,'montage');title(['  image        ','prediction']);
end
    


% synthetic nano-particle cluster image generation
% Author: Wang Yunfeng
% Date: 2020/05/26
% Southwest university 
parameter = [100;890;1280;44;2;10;3.5;87];
% nCircles m=890 n=1280 psize varity gap light gapcolor
datagenerator

clc,close all;
clearvars -except img objnum parameter
% load original img
% img = imread('data/sample1.jpg');
% remove the bottom black bar
img1 = img;
if ndims(img1) == 3
    img2 = im2double(rgb2gray(img1));
else
    img2 = im2double(img1);
end
[m,n]=size(img2);

%%
% select target particle using a bounding box
figure,imshow(img2);
set(gcf,'outerposition',get(0,'screensize'));
target_rect = imrect;
target_pos = round(getPosition(target_rect));
% crop from original image
target = imcrop(img2,target_pos);
delete(target_rect);close(figure(gcf));
figure,imshow(target);
set(gcf,'outerposition',get(0,'screensize'));
% draw a poly to cover the target particle region
roi_poly = impoly;
roi_mask = createMask(roi_poly);
seg = double(roi_mask).*target;
seg_props = regionprops(double(roi_mask),{'Centroid'});
pos_center = seg_props.Centroid;
pos_val = getPosition(roi_poly);
pos_r = round(pos_val);
pos = zeros(size(seg));
for i = 1:length(pos_r)
    pos(pos_r(i,2),pos_r(i,1))=256-i;
end
[posm,posn] = find(pos==max(pos(:)));
pos([posm-1:posm+1],[posn-1:posn+1])=255;

% kernel = [0 0.2 0;0.2 0.8 0.2;0 0.2 0];
% pos = im2double(conv2(pos,kernel,'same'));
delete(roi_poly);
close(figure(gcf));
clearvars target_rect roi_poly

ptc = seg;
ptc_mask = double(roi_mask);
figure,imshowpair(ptc,double(ptc_mask),'montage');

%%
% target particle as a seed

% adding alpha channel
% alpha = 0.1;
% props = regionprops(im2double(ptc_mask),'PixelList');
% ptc_alpha = zeros(size(ptc));
% ptc_alpha(find(ptc_mask==1)) = 1;
% se1 = strel('disk',1);
% ptc_alpha = imdilate(ptc_alpha,se1);
% ptc_alpha(find((ptc_alpha-ptc_mask)==1)) = alpha;
% figure,imshowpair(ptc_mask,ptc_alpha);
% ptc = cat(3,ptc,ptc,ptc);
% imwrite(ptc, 'ptc1.png', 'png', 'Alpha',ptc_alpha);

%save the seed
imwrite(ptc, 'data/ptc.png');
imwrite(ptc_mask,'data/ptc_mask.png');
imwrite(pos,'data/pos.png');

% imwrite(ptc,'ptc.png');

%%

clearvars -except ptc ptc_mask pos pos_val pos_center parameter
ptc = im2uint8(ptc);
% clc,clear all,close all;
% [ptc,~,ptc_alpha] = imread('ptc.png');
% ptc = double(rgb2gray(ptc));
% ptc_mask = imread('ptc_mask.png');
% figure,imshowpair(ptc,double(ptc_mask),'montage');
% imwrite(ptc_mask,'ptc_mask.png');
% imwrite(ptc,'ptc.png');
props = regionprops(im2double(ptc_mask),'centroid');
b = bwboundaries(ptc_mask);b = b{1};
centre  = props.Centroid;
dists = sqrt((b(:,1)-centre(2)).^2 + (b(:,2)-centre(1)).^2);
maxdist = max(dists);

%%
% nCircles = 800;
% m = 506;n = 728;
% psize = 16; varity = 2;
% gap = 10;
% light = 3.5;
% gapcolor = 87;

nCircles = parameter(1);
m = parameter(2);n = parameter(3);
psize = parameter(4); varity = parameter(5);
gap = parameter(6);
light = parameter(7);
gapcolor = parameter(8);
blankarea = ones(m,n);
tic;
circles = findcircleplace(nCircles,psize,maxdist,varity,gap,blankarea,light);
toc;
%%
% col 1: x-coordinate of the center position
% col 2: y-coordinate of the center position
% col 3: radius
% col 4: the scale factor
% col 5: rotatation fator
% col 6: flip horizontal or flip vertical
% col 7: light variance
%
% generate the synthetic image
% set a bunch of backgrounds
imbg = zeros(m+100,n+100);
imbg_mask = uint16(zeros(m+100,n+100));
temp = cell(nCircles,1);
pos_temp = cell(nCircles,1);

[map,map_mask,temp,pos_temp] = synimage(circles,imbg,imbg_mask,ptc,ptc_mask,pos,temp,pos_temp);
map = imcrop(map,[51,51,n-1,m-1]);
map_mask = imcrop(map_mask,[51,51,n-1,m-1]);
% fill the gap
map = uint8(map);
map(find(map==0))=gapcolor;


%%

% imshow
colormap = label2rgb(map_mask);
figure;imshowpair(map,colormap,'montage');
figure;imshowpair(map,colormap);
% map(find(map == 0))=100;
figure,imshow(map);

% mask_props = regionprops('table',map_mask,{'Centroid','BoundingBox','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength'});
% BoundingBox = mask_props.BoundingBox;
% Eccentricity = mask_props.Eccentricity;
save('data/circle.mat','map','map_mask','temp');
save('data/pos.mat','pos_temp','pos_val','pos_center');

% colormap = label2rgb(map_cir);panbie
% figure;imshowpair(colormap,map_cir,'montage');

%%

function [map,map_mask,temp,pos_temp]= synimage(circles,imbg,imbg_mask,ptc,ptc_mask,pos,temp,pos_temp)
[m,n] = size(imbg);
%make a synthesized image
circles(:,1) = circles(:,1)+50;
circles(:,2) = circles(:,2)+50;
% map_cir = uint16(zeros(m+100,n+100));
map_mask = imbg_mask;
map = imbg;
%[mm, nn] = meshgrid(1:n+100, 1:m+100);
for i = 1:length(circles)
    %map_cir((nn - circles(i,1)).^2 + (mm - circles(i,2)).^2 <= circles(i,3).^2) = uint16(i);
    ptcl = ptc;
    ptcl(find(ptc~=0)) = ptc(find(ptc~=0))-circles(i,7);
    ptc = ptcl;
    % image scale
    i;
    ptc_temp = imresize(ptc,circles(i,4),'nearest');
    ptc_mask_temp = imresize(ptc_mask,circles(i,4),'nearest');
    ptc_pos_temp = imresize(pos,circles(i,4),'nearest');
    % image rotate
    ptc_temp = imrotate(ptc_temp,circles(i,5),'crop');
    ptc_mask_temp = imrotate(ptc_mask_temp,circles(i,5),'crop');
    ptc_pos_temp = imrotate(ptc_pos_temp,circles(i,5),'crop');
    % whether fliped or not
    if circles(i,6)
        ptc_temp = flip(ptc_temp,circles(i,6));
        ptc_mask_temp = flip(ptc_mask_temp,circles(i,6));
        ptc_pos_temp = flip(ptc_pos_temp,circles(i,6));
    end
    
    % generate a cluster image
    [m_temp,n_temp] = size(ptc_temp);
    c_temp = round(circles(i,:));
    ltopx = c_temp(1)-round(n_temp/2);
    ltopy = c_temp(2)-round(m_temp/2);
    % New blank layer
    map_temp = zeros(m,n);
    map_mask_temp = uint16(zeros(m,n));
    map_pos_temp = double(zeros(m,n));
    % map layer
    map_temp([ltopx:ltopx+m_temp-1],[ltopy:ltopy+n_temp-1]) = ptc_temp;
    map(find(map_temp~=0)) = map_temp(find(map_temp~=0));
    map(find(map_temp~=0)) = map(find(map_temp~=0));
    % map mask layer
    map_mask_temp([ltopx:ltopx+m_temp-1],[ltopy:ltopy+n_temp-1]) = ptc_mask_temp;
    temp{i} = map_mask_temp([51:m-50],[51:n-50]);
    map_mask(find(map_temp~=0)) = uint16(i);
    % poly position layer
    map_pos_temp([ltopx:ltopx+m_temp-1],[ltopy:ltopy+n_temp-1]) = ptc_pos_temp;
    pos_temp{i} = im2uint8(map_pos_temp([51:m-50],[51:n-50])/255);
end
end

%%
function circles = findcircleplace(nCircles,psize,maxdist,varity,gap,blankarea,light)
% generate the properties of partical cluster
circles = zeros(nCircles,3);
[mm,nn] = size(blankarea);
% for better view
figure,hold on
set(gca,'XTick',[0:100:1300]);
set(gca,'YTick',[0:100:900]);
title(['particles placement'])
for i=1:nCircles
    %Flag which holds true whenever a new circle was found
    newCircleFound = false;
    isblank = false;
    %loop iteration which runs until finding a circle which doesnt intersect with previous ones
    while ~newCircleFound
        
        x = randi([1,nn],1);
        y = randi([1,mm],1);
        %r = 10+10*rand(1);
        r = normrnd(psize,varity);%15 2
        
        %calculates distances from previous drawn circles
        prevCirclesY = circles(1:i-1,1);
        prevCirclesX = circles(1:i-1,2);
        prevCirclesR = circles(1:i-1,3);
        distFromPrevCircles = ((prevCirclesX-x).^2+(prevCirclesY-y).^2).^0.5;
        
        %if the distance is not to small - adds the new circle to the list
        if i==1 || sum(distFromPrevCircles<=(r+prevCirclesR-gap))==0
            newCircleFound = true;
            circles(i,[1,2,3]) = [y x r];
            circle(x,y,r)
        end
    end
end
hold off;
% set factor controlling the particle
circles(:,4) = circles(:,3)/maxdist;
circles(:,5) = randi([-90 90],nCircles,1);
circles(:,6) = randi([0 2],nCircles,1);
circles(:,7) = sort(0.03*normrnd(0,light,nCircles,1),'descend');

end

%%
function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
end














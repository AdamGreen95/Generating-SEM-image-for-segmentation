% synthetic nano-particle cluster image generation
% Author: Wang Yunfeng
% Date: 2020/07/05
% Southwest university

cd /home/hadmatter/桌面/homogeneous_nanoparitcle_cluster
close all,clc,clear;
%load all files
dvar = matfile('Datasets/data.mat');

imnum = 1;% 1-200
img = reshape(dvar.img(imnum,:),890,1280);
ptcnum = ptcplot(img,dvar,imnum); %show all the particle imgs
ptcnum = 2;

%%

ptc = cell2mat(dvar.ptcell(imnum,2*ptcnum-1));
ptc_mask = cell2mat(dvar.ptcell(imnum,2*ptcnum));

global m n maxdist kernel_s
m = 890; n = 1280; % image size
kernel_s = 21; % >=21

parameter1 = [200;3;15;1];
parameter2 = [parameter1(1:2);parameter1(4)];

shading = 1; % select 0 or 1


%% placement rules
props = regionprops(ptc_mask,'centroid');
b = bwboundaries(ptc_mask);b = b{1};
centre  = props.Centroid;
dists = sqrt((b(:,1)-centre(2)).^2 + (b(:,2)-centre(1)).^2);
maxdist = max(dists);% average radius of the particle


tic
% The coincidence of any two circles cannot exceed 20 pixels
circles1 = findcircleplace1(parameter1);
toc
%
tic
% place the circles randomly
circles2 = findcircleplace2(parameter2);
toc
% circles1 & circles2:
% col 1: x-coordinate of the center position
% col 2: y-coordinate of the center position
% col 3: radius
% col 4: the scale factor
% col 5: rotatation fator
% col 6: flip horizontal or flip vertical
% col 7: light variance

%% generate the synthetic image

% set a bunch of backgrounds
% backgrounds padding with 200 pixels
% 100 pixels for up, down, left, and right  boundariies

imbg = zeros(m+400,n+400);
se = strel('disk',100);
imgbr = imfilter(img,fspecial('gaussian', [kernel_s kernel_s], kernel_s),'same');
imgo = 0.5*imopen(im2double(img),se)+0.5*im2double(imgbr);
imbg(201:end-200,201:end-200) = imgo;

imbg = zeros(m+400,n+400);

imbg_mask = uint16(zeros(m+400,n+400));

tic
[map,map_mask] = synimage(circles1,imbg,imbg_mask,ptc,ptc_mask,shading);
map2 = imadjust(map,[],[min(im2double(img(:))),1]);
map3 = imfilter(map2,fspecial('gaussian',[3 3],3),'same');

[map_random,map_mask_random] = synimage(circles2,imbg,imbg_mask,ptc,ptc_mask,0);
%[map1,map_mask1] = synimage(circles2,imbg,imbg_mask,ptc,ptc_mask,shading);
toc
%%
figure;
subplot(2,2,[1,3]),imshow(img)
subplot(2,2,2),imshow(map3)
subplot(2,2,4),imshow(map_random)
% subplot(2,2,2),imshowpair(map,map_mask,'montage')
% subplot(2,2,4),imshowpair(map1,map_mask1,'montage')
parameter = parameter1(1:3);

%% function file

function ptcnum = ptcplot(img,dvar,num)
% show all the particle imgs(2~5) and its mask in one image
pvarind = cellfun(@isempty,dvar.ptcell(num,:));
plotnum = length(find(pvarind==0));
ptcnum = randi([1,plotnum/2],1);
figure;
set(gcf,'outerposition',get(0,'screensize'));
for i = 1:plotnum
    if mod(i,2)==0
        j = (i+plotnum)/2;
        subplot(5,plotnum/2,j)
        imshow(cell2mat(dvar.ptcell(num,i)))
    else
        j = (i+1)/2;
        subplot(5,plotnum/2,j)
        imshow(cell2mat(dvar.ptcell(num,i)))
    end
    subplot(5,plotnum/2,[plotnum+1:plotnum*5/2])
    imshow(img)
end
end

function circles = findcircleplace1(parameter)
% generate the properties of partical cluster
global m n maxdist
psize = maxdist;
blankarea = ones(m,n);

nCircles = parameter(1);
varity = parameter(2);
gap = parameter(3);
light = parameter(4);

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
    i
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

function circles = findcircleplace2(parameter);

% generate the properties of partical cluster
global m n maxdist
psize = maxdist;
nCircles = parameter(1);
varity = parameter(2);
light = parameter(3);
circles = zeros(nCircles,3);

figure,hold on
set(gca,'XTick',[0:100:1300]);
set(gca,'YTick',[0:100:900]);
title(['particles placement'])
for i=1:nCircles
    x = randi([1,n],1);
    y = randi([1,m],1);
    r = normrnd(psize,varity);
    circles(i,[1,2,3]) = [y x r];
    circle(x,y,r)
end
hold off;
circles(:,4) = circles(:,3)/maxdist;
circles(:,5) = randi([-90 90],nCircles,1);
circles(:,6) = randi([0 2],nCircles,1);
circles(:,7) = sort(0.03*normrnd(0,light,nCircles,1),'descend');
end

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

function [map,map_mask]= synimage(circles,imbg,imbg_mask,ptc,ptc_mask,shading)
[m,n] = size(imbg);
%make a synthesized image
circles(:,1) = circles(:,1)+200;
circles(:,2) = circles(:,2)+200;
map_mask = imbg_mask;
map = imbg;
for i = 1:length(circles)
    ptcl = ptc;
    ptcl(find(ptc~=0)) = ptc(find(ptc~=0))-0.05*circles(i,7);% light change
    ptc = ptcl;
    % image scale
    i;
    ptc_temp = imresize(ptc,circles(i,4),'nearest');
    ptc_mask_temp = imresize(ptc_mask,circles(i,4),'nearest');
    % image rotate
    ptc_temp = imrotate(ptc_temp,circles(i,5),'crop');
    ptc_mask_temp = imrotate(ptc_mask_temp,circles(i,5),'crop');
    % whether fliped or not
    if circles(i,6)
        ptc_temp = flip(ptc_temp,circles(i,6));
        ptc_mask_temp = flip(ptc_mask_temp,circles(i,6));
    end
    % generate a cluster image
    [m_temp,n_temp] = size(ptc_temp);
    c_temp = round(circles(i,:));
    ltopx = c_temp(1)-round(n_temp/2);
    ltopy = c_temp(2)-round(m_temp/2);
    % New blank layer
    map_temp = zeros(m,n);
    map_mask_temp = uint16(zeros(m,n));
    % map layer
    map_temp([ltopx:ltopx+m_temp-1],[ltopy:ltopy+n_temp-1]) = ptc_temp;
    map(find(map_temp~=0)) = map_temp(find(map_temp~=0));
    
    % map mask layer
    map_mask_temp([ltopx:ltopx+m_temp-1],[ltopy:ltopy+n_temp-1]) = ptc_mask_temp;
    map_mask(find(map_temp~=0)) = uint16(i);
    
    % shading layer
    if shading==1
        shadL_temp = shadinglayer(map_mask_temp);
        shadL_temp(find(map_mask_temp~=0)) = 0;
        map = map-0.5*double(shadL_temp);
    end
    
end
map = imcrop(map,[201,201,n-400,m-400]);
map_mask = imcrop(map_mask,[201,201,n-400,m-400]);
end

function sd_mask = shadinglayer(mask)
global kernel_s
mask = gpuArray(mask);
sh_b1 = double(edge(mask));
sh_b2 = imfilter(sh_b1,fspecial('gaussian',[9 9],9),'same');
sh_b = sh_b2;
sh_b(find(sh_b2~=0))=1;
Ig1 = imfilter(sh_b,fspecial('gaussian',[9 9],9),'same');
Ig2 = imfilter(sh_b,fspecial('gaussian',[kernel_s-10 kernel_s-10],kernel_s-10),'same');
Ig3 = imfilter(sh_b,fspecial('gaussian',[kernel_s kernel_s],kernel_s),'same');
Ig4 = imfilter(sh_b,fspecial('gaussian',[3*kernel_s+10 3*kernel_s+10],3*kernel_s+10),'same');
Ig5 = imfilter(sh_b,fspecial('gaussian',[101 101],101),'same');
co = [max(Ig5(:)),max(Ig4(:)),max(Ig3(:)),max(Ig2(:)),max(Ig1(:))];
co = co/sum(co);
Ig = co(1)*Ig1+co(2)*Ig2+co(3)*Ig3+co(4)*Ig4+co(5)*Ig5;
sd_mask = gather(Ig);
end








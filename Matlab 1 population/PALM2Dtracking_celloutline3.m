function PALM2Dtracking_celloutline3(folder)
close
folder_name=folder;
%% find the location of the bright field movie
% folder_name='E:/tamu/tamufit/PALM Nan lab/FrzCD/08142017FrzCDGFPcephalexin8htreatment7';
bf_folder='bfbefore';
d3=dir([folder_name '/' bf_folder]);
for i=1:1:length(d3)
    if strcmp(d3(i).name,'bfbefore_MMStack_Pos0.ome.tif')==1
        bf_movie=d3(i).name;
    end
end

img=65535-double(imread([folder_name '/' bf_folder '/' bf_movie]));
% img=double(imread([folder_name '/' bf_folder '/' bf_movie]));
cell_num=load([folder_name '/' bf_folder '/'  'mcell_num.txt']);
for n=1:1:cell_num
%% roi for the spot
xy=load([folder_name '/' bf_folder '/mask' num2str(n) '.txt']);
row_low=round(min(xy(:,2)));
col_low=round(min(xy(:,1)));
row_high=round(max(xy(:,2)));
col_high=round(max(xy(:,1)));
img2=img(row_low:row_high,col_low:col_high);
try
%set(gcf,'position',get(0,'screensize'));
subplot(1,2,1)
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
title([num2str(n) ' / ' num2str(cell_num)])
%% threshold
img3=mat2gray(img2);
level = graythresh(img3);
bw = im2bw(img3,level+0.15);
subplot(1,2,2)
imshow(bw);
bw2=bwselect(bw,8);
bw3=edge(bw2,'log');
%% find the edge
% bw2=edge(img2,'log');
% % bw2=edge(img,'sobel');
% % bw2=edge(img,'roberts');
% % bw2=edge(img,'prewitt');
% subplot(1,2,2)
% imshow(bw2)
% pause(1)
% close
%% connect the outline
% bw2=double(bw2);
% set(gcf,'position',get(0,'screensize'));
% imshow(bw2,'DisplayRange',[min(min(bw2)),max(max(bw2))],'InitialMagnification','fit')
% button=1;
% while button==1
% [xp,yp,bt]=ginput(1);
% xp=round(xp);
% yp=round(yp);
% button=bt;
% if button==1
% bw2(yp,xp)=1;
% imshow(bw2,'DisplayRange',[min(min(bw2)),max(max(bw2))],'InitialMagnification','fit')
% title('connect the outline')
% end
% end
%% select the outline
% set(gcf,'position',get(0,'screensize'));
% imshow(bw2,'DisplayRange',[min(min(bw2)),max(max(bw2))],'InitialMagnification','fit')
% title('select the outline')
% [p,xi,yi]=roipoly;
% bw3=bw2.*p;
% %% dilate the outline
se=strel('square',3);
bw4=imdilate(bw3,se);
%% obtain the coordinates of outline
b=bwboundaries(bw4);
%% check the outline
set(gcf,'position',get(0,'screensize'));
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
hold on
plot(b{1}(:,2),b{1}(:,1),'r')
xy_outline=[b{1}(:,2),b{1}(:,1)];
save([folder_name '/' bf_folder '/xy_outline' num2str(n) '.txt'],'-ASCII','-TABS','xy_outline');
pause(1)
close
catch me
set(gcf,'position',get(0,'screensize'));
subplot(1,2,1)
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
%% threshold
img3=mat2gray(img2);
level = graythresh(img3);
bw = im2bw(img3,level+0.15);
subplot(1,2,2)
imshow(bw);
bw2=bwselect(bw,8);
bw3=edge(bw2,'log');    
se=strel('square',3);
bw4=imdilate(bw3,se)
%% obtain the coordinates of outline
b=bwboundaries(bw4);
%% check the outline
set(gcf,'position',get(0,'screensize'));
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
hold on
plot(b{1}(:,2),b{1}(:,1),'r')
xy_outline=[b{1}(:,2),b{1}(:,1)];
save([folder_name '/' bf_folder '/xy_outline' num2str(n) '.txt'],'-ASCII','-TABS','xy_outline');
pause(1)
close
end
end
%% thresholding
% img3=mat2gray(img2);
% level = graythresh(img3);
% bw = im2bw(img3,level);
% imshow(bw);
end
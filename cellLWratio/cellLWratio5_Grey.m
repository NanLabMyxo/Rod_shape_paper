clc
clear
folder_name='/Users/chen/Desktop/cellLWratio/20180723/mrec/dz2 spores';
key_word='dz2 spores';
d3=dir(folder_name);
for i=3:1:length(d3)
d3(i).name;
if contains(d3(i).name,key_word)==1
sub_folder=d3(i).name;
bf_folder=[folder_name,'/' sub_folder];
d4=dir(bf_folder);
cell_ratio=[];
w=0;
for ii=3:1:length(d4)
bf_file=d4(ii).name;
if contains(bf_file,'data')==1
    break;
else
img=imread([bf_folder '/' bf_file]);
img2=255-rgb2gray(img);
%%
img3=mat2gray(img2);
level = graythresh(img3);
bw2 = im2bw(img3,level);
%% dilate the outline
se=strel('square',3);
bw=imdilate(bw2,se);
%%
% subplot(1,2,1)
% imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
% set(gcf,'position',get(0,'screensize'));
% title(bf_folder)
% subplot(1,2,2)
% imshow(bw)
% pause(2)
% close
%%
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
set(gcf,'position',get(0,'screensize'));
title(bf_folder)
% pause(5)
%%
b=bwboundaries(bw,'noholes');
for j=1:1:length(b)
b2=b{j};
bx=b2(:,2);
by=b2(:,1);
if length(b2)>200 && length(b2)<800  %% lenth range and repeat number
hold on
plot(bx,by)
pause(1)
% length(b2)
%% calculate the length to width ratio
w=w+1;
bws=bw(min(by):max(by),min(bx):max(bx));
region_length=regionprops(bws,'MajorAxisLength');
region_width=regionprops(bws,'MinorAxisLength');
cell_long=[];
% cell_short=[];
for v=1:1:length(region_length)
    cell_long(v)=region_length(v).MajorAxisLength;
%     cell_short(v)=region_width(v).MinorAxisLength;
end
[u,v]=max(cell_long);
cell_length=u;
cell_width=region_width(v).MinorAxisLength;
cell_ratio(w)=cell_length/cell_width;
end
end %% for j+1:1:length(b)
pause(2)
hold off
end
end  %% d4
cell_ratio=cell_ratio';
save([bf_folder '/data.txt'],'-ASCII','-TABS','cell_ratio')
end  %% for each folder
end  %% all folders
pause(2)
close

%% for folder contains multiple images
clc
clear
folder_name='S:\Nan_Lab\Autumn\2021-03\20210317 sporulation d4628\20210317 sporulation d4628 2h';
cd(folder_name)
key_word='2021';
d3=dir('*.tif');
for i=1:1:length(d3)
d3(i).name;
if ~isempty(strfind(d3(i).name,key_word)) && ~isempty(strfind(d3(i).name,'tif'))
img2=65535-double(imread([folder_name '/' d3(i).name]));
img=img2; %%
bw=edge(img2,'log');
%%
% subplot(1,2,1)
% imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
% set(gcf,'position',get(0,'screensize'));
% subplot(1,2,2)
% imshow(bw)
% pause(1)
% close
%%
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
set(gcf,'position',get(0,'screensize'));
b=bwboundaries(bw,'noholes');
title ([num2str(i) '/' num2str(length(d3))])
%% sort the movie and save
s1=d3(i).name(1:end-7)
if i==1
s0=s1;
w=0;  
cell_ratio=[];
cell_length=[];
cell_width=[];
end
if strcmp(s0,s1)
else
cell_ratio=cell_ratio';
cell_length=cell_length';
cell_width=cell_width';
save([folder_name '/' s0 '_ratio.txt'],'-ASCII','-TABS','cell_ratio')
save([folder_name '/' s0 '_length.txt'],'-ASCII','-TABS','cell_length')
save([folder_name '/' s0 '_width.txt'],'-ASCII','-TABS','cell_width')
s0=s1;
w=0;  
cell_ratio=[];
cell_length=[];
cell_width=[];
end
for j=1:1:length(b)
b2=b{j};
bx=b2(:,2);
by=b2(:,1);
if length(b2)>30 && length(b2)<90  && by(3)~=by(end-2) %% lenth range and repeat number
[m,n]=max(bx);
bx1=bx(1:n);
bx2=bx(n:end);
bx1d=diff(bx1);
bx2d=diff(bx2);
if min(bx1d)>=0 && max(bx2d)<=0 %% x number in order
n2=find(bx==max(bx));
n3=find(by==max(by));
if  n2(1)-1>0 && n3(1)-1>0 && by(n2(1)-1)~=by(n2(end)+1)&& bx(n3(1)-1)~=bx(n3(end)+1)  %% repeat number in around maximum of y
hold on
plot(bx,by)
% pause(1)
% close
%% calculate the length to width ratio
w=w+1;
bws=bw(min(by):max(by),min(bx):max(bx));
% imshow(bws)
% set(gcf,'position',get(0,'screensize'));
% pause(1)

region_length=regionprops(bws,'MajorAxisLength');
region_width=regionprops(bws,'MinorAxisLength');
cell_long=[];
% cell_short=[];
for v=1:1:length(region_length)
    cell_long(v)=region_length(v).MajorAxisLength;
%     cell_short(v)=region_width(v).MinorAxisLength;
end
[u,v]=max(cell_long);
cell_length(w)=u;
cell_width(w)=region_width(v).MinorAxisLength;
cell_ratio(w)=cell_length(w)/cell_width(w);
end
end
end
end %% for j+1:1:length(b)
end  %% for each image
pause(0.5)
end  %% all images
if i==length(d3)
cell_ratio=cell_ratio';
cell_length=cell_length';
cell_width=cell_width';
save([folder_name '/' s1 '_ratio.txt'],'-ASCII','-TABS','cell_ratio')
save([folder_name '/' s1 '_length.txt'],'-ASCII','-TABS','cell_length')
save([folder_name '/' s1 '_width.txt'],'-ASCII','-TABS','cell_width')
end

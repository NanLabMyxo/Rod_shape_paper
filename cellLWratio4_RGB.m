clc
clear
folder_name='/Users/chen/Desktop/cellLWratio/20180723/mrec/dz2 spores';
key_word='dz2';
d3=dir(folder_name);
for i=3:1:length(d3)
d3(i).name;
if ~isempty(strfind(d3(i).name,key_word)) 
%if isempty(strfind(d3(i).name,'Stack'))
sub_folder=d3(i).name;
bf_folder=[folder_name,'/' sub_folder];
bf_file=[sub_folder '_MMStack_Pos0.ome'];
img2=65535-double(imread([bf_folder '/' bf_file '.tif']));
img=img2; %%
bw=edge(img2,'log');
%%
% subplot(1,2,1)
% imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
% set(gcf,'position',get(0,'screensize'));
% title(bf_folder)
% subplot(1,2,2)
% imshow(bw)
% pause(1)
% close
%%
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
set(gcf,'position',get(0,'screensize'));
title(bf_folder)
b=bwboundaries(bw,'noholes');
w=0;
cell_ratio=[];
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
cell_length=u;
cell_width=region_width(v).MinorAxisLength;
cell_ratio(w)=cell_length/cell_width;
end
end
end
end %% for j+1:1:length(b)
pause(1)
hold off
cell_ratio=cell_ratio';
save([bf_folder '/data.txt'],'-ASCII','-TABS','cell_ratio')


end  %% for each folder
end  %% all folders
pause(2)
close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% combine data
k=1;
directory=struct([]);
for i=3:1:length(d3)
    if isempty(strfind(d3(i).name,key_word))
    else
        sub_folder=d3(i).name;
directory(k).name=sub_folder(1:end);
k=k+1;
    end
end
%%
%% 
first=1;
d=directory(1).name; %% first folder
b=regexp(d(end-1:end),'\d'); %% last two digit number
n=length(d)-length(b);
if length(directory)==1
last=1;
sub_folder=directory(1).name;
folder=[folder_name,'/' sub_folder];
cellLWratio4combine(folder,first,last)
end
w=1;
for i=1:1:length(directory)-1
if strncmp(directory(i).name,directory(i+1).name,n)==1
w=w+1;
if i==length(directory)-1
last=w
sub_folder=directory(i).name;
folder=[folder_name,'/' sub_folder];
cellLWratio4combine(folder,first,last)
end
else
last=w
sub_folder=directory(i).name;
folder=[folder_name,'/' sub_folder]
pause(3)
cellLWratio4combine(folder,first,last)
w=1;
d=directory(i+1).name;
b=regexp(d(end-1:end),'\d');
n=length(d)-length(b);
end %% strncmp

end %% for i=1:1:length(directory)

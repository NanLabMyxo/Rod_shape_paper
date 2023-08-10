function PALM2Dtracking_cellselection4(folder)
close
%% total cell number
folder_name=folder;
%% find the location of the bright field movie
bf_folder='bfbefore';
d3=dir([folder_name '/' bf_folder]);
for i=1:1:length(d3)
    if strcmp(d3(i).name,'bfbefore_MMStack_Pos0.ome.tif')==1
        bf_movie=d3(i).name;
    end
end
%% read the bright field movie
img=double(imread([folder_name '/' bf_folder '/' bf_movie]));
%set(gcf,'position',get(0,'screensize'));
imshow(img,'DisplayRange',[min(min(img)),max(max(img))],'InitialMagnification','fit')
title([folder_name '/' bf_folder])
%% select ROI regions
p1=zeros(size(img));
button=1;
i=0;
while button==1
keydown = waitforbuttonpress;
if (keydown == 1)
break;
end
i=i+1;
[p,xi,yi]=roipoly;
p1=p1+p;
[k,w]=find(p==1);
hold on
text(w(1),k(1),num2str(i),'color','r','fontsize',20);
% imagesc(p)
pause(1)
% close
p=uint16(p);
imwrite(p,[folder_name '/' bf_folder '/mask' num2str(i) '.tif'],'tiff','WriteMode','overwrite' )
xy=[xi,yi];
save([folder_name '/' bf_folder '/mask' num2str(i) '.txt'],'-ascii','-TABS','xy');
end
hold off
imagesc(p1)
pause(2)
close
p2=img.*p1;
set(gcf,'position',get(0,'screensize'));
imshow(p2,'DisplayRange',[min(min(p2)),max(max(p2))],'InitialMagnification','fit')
pause(2)
close
%% save the mask
p1=uint16(p1);
imwrite(p1,[folder_name '/' bf_folder '/mask.tif'],'tiff','WriteMode','overwrite' )
save([folder_name '/' bf_folder '/'  'mcell_num.txt'],'-ascii','-TABS','i');
end


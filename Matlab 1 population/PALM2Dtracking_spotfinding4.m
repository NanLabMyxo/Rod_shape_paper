function PALM2Dtracking_spotfinding4(folder,first_frame,last_frame,pixel_size,values)
close
%% parameters
movie_start=first_frame;
movie_end=last_frame;
pixel_size_standard=0.16; %% um
r=fix(2*pixel_size_standard/pixel_size); %% spot size 2*r +1 by 2*r+1
spot_num=5;
folder_name=folder;
value=values;
%% find the location of bright field mask 
% folder_name='E:/tamu/tamufit/PALM Nan lab/Mray/20170620mraypamcherrya22treatment6';
d=dir(folder_name);
bb=[];
k=1;
for i=3:1:length(d)
    if isempty(strfind(d(i).name,'bf')) && isempty(strfind(d(i).name,'DS_Store'))
        bb(k)=d(i).isdir+0;
        k=k+1;
    end
end
bf_folder='bfbefore';
%% find the location of fluorescent movie
for i=3:1:length(d)
    if isempty(strfind(d(i).name,'bf')) && isempty(strfind(d(i).name,'DS_Store'))
        fl_folder=d(i).name;
        break;
    end
end
%%
cell_num=load([folder_name '/' bf_folder '/'  'mcell_num.txt']);
for n=1:1:cell_num
%% load the mask
xy=load([folder_name '/' bf_folder '/mask' num2str(n) '.txt']);
%% roi for the spot
row_low=round(min(xy(:,2)));
col_low=round(min(xy(:,1)));
row_high=round(max(xy(:,2)));
col_high=round(max(xy(:,1)));
%% check the roi region
% imshow(bf_mask,'DisplayRange',[min(min(bf_mask)),max(max(bf_mask))],'InitialMagnification','fit')
% hold on
% line(xy(:,1),xy(:,2),'color','r')
% pause(2)
%% find the spots for all the movies
for w=1:1:length(bb)
if value==0
h=waitbar((w+(n-1)*length(bb))/(cell_num*length(bb)),' please wait...');
end
fl_movie=[fl_folder(1:end-1) num2str(w) '_MMStack_Pos0.ome.tif'];
%% find the spots
x_sum=[];
y_sum=[];
frame_sum=[];
for i=movie_start:1:movie_end
img=double(imread([folder_name '/' fl_folder(1:end-1) num2str(w) '/' fl_movie],i));
img2=img(row_low:row_high,col_low:col_high);
img3=img(row_low:row_high,col_low:col_high); %% for check
peak_low=1.1*mean(mean(img2));
ave_low=1.1*mean(mean(img2));
for ii=1:1:spot_num
[pss1,pii]=max(img2);
[pss2,piii]=max(pss1); %% pss2 is the maximal intensity
yp=pii(piii);
xp=piii;
[sm,sn]=size(img2);
if yp-r<1 || xp-r<1 || yp+r>sm || xp+r>sn
    break
else
spot_region=double(img2(yp-r:yp+r,xp-r:xp+r));%% for threshold
img2(yp-r:yp+r,xp-r:xp+r)=zeros(2*r+1,2*r+1); %% remove found spots
end
%  max(max(spot_region))
if max(max(spot_region))> peak_low && mean(mean(spot_region))> ave_low  
%% Gaussian fitting
xp;
yp;
pm=max(max(spot_region));
bm=mean(mean(img2));
g0=[bm pm-bm r+1 r+1 1];
options = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolX',1e-10,'LargeScale','on');
[x1,resnorm,residual,exitflag,output,lambda] = lsqcurvefit(@gau,double(g0),r,double(spot_region),[], [],options); 
x1(4); % x coordinate
x1(5); % y coordinate
x1(3);%% size of spots
x0=xp-r+x1(4)-1;
y0=yp-r+x1(5)-1;
% if sqrt((x0-xp)^2+(y0-yp)^2)<1
% if x1(4)>1 && x1(5)>1 && x1(4)<4 && x1(5)<4 && x1(3)<3
%% check the spots found
if value==1
set(gcf,'position',get(0,'screensize'));
imshow(img3,'DisplayRange',[min(min(img3)),max(max(img3))],'InitialMagnification','fit')
title([ 'moive' num2str(w) '   ' num2str(i) '/' num2str(movie_end)])
hold on
plot(xp, yp,'b*');
plot(x0, y0,'r*');
% imagesc(spot_region2)
pause(1)
end
%% combine all data
x_sum=[x_sum,x0];
y_sum=[y_sum,y0];
frame_sum=[frame_sum,i];
% end
end   
end
end
%% save the data
xy_sum=[x_sum',y_sum',frame_sum'];
save([folder_name '/' fl_folder(1:end-1) num2str(w) '/mask ' num2str(n) ' xy_sum' num2str(w) '.txt'],'-ASCII','-TABS','xy_sum');
if value==0
close(h)
end
end
close
end
function f=gau(g,r)
[x,y]=meshgrid(1:1:2*r+1);
f=g(1)+g(2)*exp((-(x-g(4)).^2-(y-g(5)).^2)/(2*g(3).^2));
end
end


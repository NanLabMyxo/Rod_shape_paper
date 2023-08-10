function PALM2Dtracking_sortingautomatically9(folder,pixel_sizes,exposure_time,trajectory_length_low,trajectory_length_high,step_long,step_short,values,first_movie,last_movie)
%% calculate the diffusion coefficient from all trajectories and bootstrap
close
%% parameters
pixel_size=pixel_sizes; %% um
exp_time=exposure_time; %% second
tj_len_low=trajectory_length_low; %% trajectory length longer than this will be used, number of positions in each trajectory (at least 4 points)
tj_len_high=trajectory_length_high;
tj_step_long=step_long; %% 2*160nm/100ms=3.2um/s
tj_step_short=step_short;
region_size=0.1; %% um
value=values;
movie_first=first_movie;
movie_last=last_movie;
MSD_matrixb=zeros(1,1);
krowb=1;
krowd=1;
%% find the location of data
folder_name=folder;
% folder_name='E:\tamu\tamufit\PALM Nan lab\FrzCD\08142017FrzCDGFPcephalexin8htreatment3';
d=dir(folder_name);
% bb=[];
% k=1;
% for i=3:1:length(d)
%     if isempty(findstr(d(i).name,'bf'))
%         bb(k)=d(i).isdir+0;
%         k=k+1;
%     end
% end
% length(bb);
bf_folder='bfbefore';
for i=3:1:length(d)
    if isempty(strfind(d(i).name,'bf')) && isempty(strfind(d(i).name,'DS_Store')) 
        fl_folder=d(i).name;
        break;
    end
end
cell_num=load([folder_name '/' bf_folder '/'  'mcell_num.txt']);
for n=1:1:cell_num
    close
if value==0
    h=waitbar(n/cell_num,'please wait...');
end
%% load the mask
xy=load([folder_name '/' bf_folder '/mask' num2str(n) '.txt']);
xy_outline=load([folder_name '/' bf_folder '/xy_outline' num2str(n) '.txt']);
%% roi for the spot
row_low=round(min(xy(:,2)));
col_low=round(min(xy(:,1)));
row_high=round(max(xy(:,2)));
col_high=round(max(xy(:,1)));
%% for outline plot
x0=xy_outline(:,1);
y0=xy_outline(:,2);
x_dif=max(x0)-min(x0);
y_dif=max(y0)-min(y0);
xm=[];
ym=[];
%% obtain the middle line
if x_dif>y_dif %% cell long axis along x direction
    k=1;
for i=min(x0):1:max(x0)
  xc=find(x0==i);    
  ym(k)= mean(y0(xc));
  xm(k)=i;
  k=k+1;
end   
% fit the middle line to polynominal function
pf1=polyfit(xm,ym,2);  
xf1=linspace(min(xm),max(xm));  
yf1=polyval(pf1,xf1);  
% plot(xm,ym,'*',xf1,yf1,'g')
end
%
if x_dif<=y_dif  %% cell long axis along y direction
    k=1;
for i=min(y0):1:max(y0)
  yc=find(y0==i) ;
  xm(k)= mean(x0(yc));
  ym(k)=i;
  k=k+1;
end   
pf1=polyfit(ym,xm,2);  
yf1=linspace(min(ym),max(ym));  
xf1=polyval(pf1,yf1);  
% plot(xm,ym,'*',xf1,yf1,'g')
end
%% movie number
for w=movie_first:1:movie_last
%% read the data
xy_sum=load([folder_name '/' fl_folder(1:end-1) num2str(w) '/mask ' num2str(n) ' xy_sum' num2str(w) '.txt']);
if isempty(xy_sum)==1
else
x_sum=xy_sum(:,1)';
y_sum=xy_sum(:,2)';
frame_sum=xy_sum(:,3)';
%% find link between two successive frame
frame_sum2=[0,frame_sum,999];
link=diff(frame_sum2);
link2=find(link>0);
link3=link2(1:end-1); %% movie frame number
tj_index=[];
tj_sum=[];
for j=1:1:length(link3)-2
for jj=link3(j):1:link3(j+1)-1 %% for frame link3(j), each time, each spot
x1=x_sum(jj);
y1=y_sum(jj);
num1=frame_sum(jj);
xy_dis=[];
frame_index=[];
k=1;
for jjj=link3(j+1):1:link3(j+1+1)-1 %% for frame link3(j+1) 
x2=x_sum(jjj);
y2=y_sum(jjj);
xy_dis(k)=sqrt((x2-x1).^2+(y2-y1).^2); %% calculate the distance for all spots in this frame
k=k+1;
end
num2=frame_sum(jjj);
[dm,dn]=min(xy_dis);
if dm<=tj_step_long && dm>tj_step_short && num2-num1==1 %%  select trajecotry based on the step distance 
tj_index=[jj,link3(j+1)+dn-1];
tj_sum=[tj_sum,tj_index]; %% frame number for trajectory in successive frame
end
end
end
%% link the trajectory
tj1=tj_sum(1:2:end);
tj1b=tj_sum(1:2:end);
tj2=tj_sum(2:2:end);
kk=1;
tj.frame={};
tj.x={};
tj.y={};
for t=1:1:length(tj1)-2
    t_sum=t;
    for p=t:1:length(tj1)
    if tj1(t)>0
    tn=find(tj1==tj1(t));
    tn2=find(tj1==tj2(tn));
    tj1(t)=0;
    if isempty(tn2)==1
    else
    t=tn2;
    t_sum=[t_sum,tn2];
    end
    end    
    end
tj_num=[tj1b(t_sum), tj2(t_sum(end))];
if length(tj_num)>tj_len_low-1 && length(tj_num)<tj_len_high+1
tj.frame(kk)={tj_num};
kk=kk+1;
end
end
%% combine all trajectory data to a structure variable tj and check and  sort
for pp=1:1:kk-1
x_tj=x_sum(tj.frame{pp});
y_tj=y_sum(tj.frame{pp});
frame_tj=frame_sum(tj.frame{pp});
tj.x(pp)={x_tj};
tj.y(pp)={y_tj};
%% remove the spots outside of the cells
x_tj0=round(x_tj);
y_tj0=round(y_tj);
x_tj1=x_tj;
y_tj1=y_tj;
s='yes';
%% remove the spots outside of the cells
if x_dif>y_dif  %% cell long axis along x direction
   if min(x_tj0)<min(x0) || max(x_tj0)>max(x0)
       s='no';
   else
   for ii=1:1:length(x_tj0)
   yc=find(x0==x_tj0(ii));
   if y_tj0(ii)< min(y0(yc)) || y_tj0(ii)> max(y0(yc))
       s='no';
   end      
   end
   end
end
%%
if x_dif<=y_dif %% cell long axis along y direction
   if min(y_tj0)<min(y0) || max(y_tj0)>max(y0)
       s='no';
   else
   for ii=1:1:length(y_tj0)
   xc=find(y0==y_tj0(ii));
   if x_tj0(ii)< min(x0(xc)) || x_tj0(ii)> max(x0(xc))
       
       s='no';
   end
   end   
   end
end
%% remove weird trajectory
% tj_steps=[];
% for i=1:1:length(x_tj)-1
% tj_steps(i)=sqrt(((x_tj(i+1)-x_tj(i))^2+(y_tj(i+1)-y_tj(i))^2));
% end
% if (max(tj_steps)-min(tj_steps))>tj_step_long
%     s='no';
% end
if strcmp(s,'yes')
%% check the trajectory on the original movie
if value==1
subplot(1,2,1)
fl_movie=[fl_folder(1:end-1) num2str(w) '_MMStack_Pos0.ome.tif'];
for jjj=1:1:length(frame_tj)
img=double(imread([folder_name '/' fl_folder(1:end-1) num2str(w) '/' fl_movie],frame_tj(jjj)));
img2=img(row_low:row_high,col_low:col_high);
set(gcf,'position',get(0,'screensize'));
% set(gcf,'position',[100 378 560 420]);
imshow(img2,'DisplayRange',[min(min(img2)),max(max(img2))],'InitialMagnification','fit')
hold on
plot(xy_outline(:,1),xy_outline(:,2),'b')  %% plot the cell contour
plot(xm,ym,'*',xf1,yf1,'g')  %% plot the middle line
plot(x_tj(1:jjj), y_tj(1:jjj),'b',x_tj(1:jjj), y_tj(1:jjj),'r*');
title(['frame ' num2str(frame_tj(jjj)) ' / movie ' num2str(w) ' / cell' num2str(n)])
pause(1)
end
%% plot
subplot(1,2,2)
plot(x_tj,y_tj,'.r',x_tj,y_tj,'b')
hold off
end
%% middle line angle
% if x_dif>y_dif 
% xcf=[min(x_tj),max(x_tj)];
% ycf=polyval(pf1,xcf);
% end
% if x_dif<=y_dif 
% ycf=[min(y_tj),max(y_tj)];
% xcf=polyval(pf1,ycf);
% end
% midline_k=(ycf(2)-ycf(1))/(xcf(2)-xcf(1));
% y-y1=k(x-x1);
% kx-y+(y1-k*x1)=0;
% a=k
% b=-1
% c=y1-k*x1
% a=midline_k;
% b=-1;
% c=ycf(1)-midline_k*xcf(1);
%% calculate the MSD
% MSD=[];
% for m=1:1:length(x_tj)-1
% MSDi=[];
% for i=1:1:length(x_tj)-m
%%%%%%%%%% calibrate the curvature of cell
% MSD_k=(y_tj(i+m)-y_tj(i))/(x_tj(i+m)-x_tj(i));
% theta=abs((MSD_k-midline_k)/(1+MSD_k*midline_k));
% MSDp=cos(atan(theta))^2*pixel_size^2*((x_tj(i+m)-x_tj(i))^2+(y_tj(i+m)-y_tj(i))^2);
% x1=x_tj(i+m);
% y1=y_tj(i+m);
% d1=pixel_size*abs(a*x1+b*y1+c)/sqrt(a^2+b^2);
% x2=x_tj(i);
% y2=y_tj(i);
% d2=pixel_size*abs(a*x2+b*y2+c)/sqrt(a^2+b^2);
% 
% a1=MSD_k;
% b1=-1;
% c1=y1-MSD_k*x1;
% xc=(c1-c)/(a-a1);
% yc=a1*xc+c1;
% if (x1-xc)*(x2-xc)>=0 %% two sports are the same side relative to the middle line
% if d1/cell_r<1 && d2/cell_r<1
% theta1=acos(d1/cell_r);
% theta2=acos(d2/cell_r);
% arc_len=abs(theta1-theta2)*cell_r;
% MSDs=arc_len^2;
% MSDi(i)=MSDp+MSDs;
% else
% MSDi(i)=pixel_size^2*((x_tj(i+m)-x_tj(i))^2+(y_tj(i+m)-y_tj(i))^2);
% end
% end
% if (x1-xc)*(x2-xc)<0 %% two sports are both sides relative to the middle line
% if d1/cell_r<1 && d2/cell_r<1
% theta1=acos(d1/cell_r);
% theta2=acos(d2/cell_r);
% arc_len=(pi-theta1-theta2)*cell_r;
% MSDs=arc_len^2;
% MSDi(i)=MSDp+MSDs;
% else
% MSDi(i)=pixel_size^2*((x_tj(i+m)-x_tj(i))^2+(y_tj(i+m)-y_tj(i))^2);
% end
% end
% end
% MSD(m)=mean(MSDi);
% end
%% no curvature calibration
MSD=[];
for m=1:1:length(x_tj)-1
MSDi=[];
for i=1:1:length(x_tj)-m
MSDi(i)=pixel_size^2*((x_tj(i+m)-x_tj(i))^2+(y_tj(i+m)-y_tj(i))^2);
end
MSD(m)=mean(MSDi);
end
%% calculate individual diffusion coefficient
if length(MSD)>3
[p1,s1]=polyfit(exp_time*(1:4),MSD(1:4),1);
% D_std= sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);
% D_std1=D_std(1)/4;
if p1(1)>0
D_coef1(krowd)=p1(1)/4;
krowd=krowd+1;
end
end

%% combine MSD 
for i=1:1:length(MSD)
MSD_matrixb(krowb,i)=MSD(i);
end
krowb=krowb+1;

end %% yes
end  %% pp
end  %% isempty
end  %% all movie
if value==0
close(h)
end
end %% all cells
close all

%% average MSD
[m,n]=size(MSD_matrixb);
for i=1:1:n
MSD_col= MSD_matrixb(:,i);
[j,k]=find(MSD_col>0);
p=1:1:length(j);
MSD_aveb(i)=mean(MSD_col(j(p)));
MSD_stdb(i)=std(MSD_col(j(p)));
end
MSD_aveb=MSD_aveb';
MSD_stdb=MSD_stdb';
D_coef1=D_coef1';
save([folder_name '/' bf_folder '/zMSD_aveb.txt'],'-ASCII','-TABS','MSD_aveb');
save([folder_name '/' bf_folder '/zMSD_stdb.txt'],'-ASCII','-TABS','MSD_stdb');
save([folder_name '/' bf_folder '/zMSD_matrixb.txt'],'-ASCII','-TABS','MSD_matrixb');
save([folder_name '/' bf_folder '/zdiffusion_coefficient_individual.txt'],'-ASCII','-TABS','D_coef1');
%% bootstrap
[m,n]=size(MSD_matrixb);
MSD_matrix1=[];
MSD_boot=[];
for w=1:1:1000
% construct a bootstrap matrix
for i=1:1:m
k=randi([1,m],1);
MSD_matrix1(i,:)=MSD_matrixb(k,:);
end
%% calculate the mean for each column
MSD_ave=[];
for i=1:1:n
MSD_col= MSD_matrix1(:,i);
[j,k]=find(MSD_col>0);
if isempty(j)
else
p=1:1:length(j);
MSD_ave(i)=mean(MSD_col(j(p)));
end
end
MSD_ave0=padarray(MSD_ave,[0 n-length(MSD_ave)],'post');
MSD_boot(w,:)=MSD_ave0;
end

MSD_bootstrap=[];
for i=1:1:n
MSD_col= MSD_boot(:,i);
[j,k]=find(MSD_col>0);
if isempty(j)
    break
else
p=1:1:length(j);
MSD_bootstrap(i)=std(MSD_col(j(p)));
end
end
MSD_bootstrap=MSD_bootstrap';
save([folder_name '/' bf_folder '/zMSD_matrixb_bootstrap.txt'],'-ASCII','-TABS','MSD_bootstrap');
%% calculate the diffusion coefficient
time_lag1=exp_time*(1:4)';
MSD1=MSD_aveb(1:4);
% g0=[1 1];
% MSD_fit=inline('g(1)+x*4*g(2)','g','x');
% options = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolX',1e-10,'LargeScale','on');
% [f,resnorm,residual,exitflag,output,lambda] = lsqcurvefit(MSD_fit,g0,time_lag1,MSD1,[],[],options); 
% f_MSD=MSD_fit(f,time_lag1);
% r_squared=1-sum((MSD1-f_MSD).^2)/sum((MSD1-mean(MSD1)).^2);
[p,s]=polyfit(time_lag1,MSD1,1);
D_std= sqrt(diag(inv(s.R)*inv(s.R')).*s.normr.^2./s.df);
D_std1=D_std(1)/4;
D_coef=p(1)/4;
D_sum=[D_coef;D_std1;m];
save([folder_name '/' bf_folder '/zdiffusion_coefficient.txt'],'-ASCII','-TABS','D_sum');
%% plot
subplot(1,2,1)
plot(exp_time*(1:1:length(MSD_aveb)),MSD_aveb,'k',exp_time*(1:1:length(MSD_aveb)),MSD_aveb,'k*');
errorbar(exp_time*(1:1:length(MSD_aveb)),MSD_aveb,MSD_bootstrap')
xlabel('Time (second)');
ylabel('MSD(\mum^{2})')
subplot(1,2,2)
bar(1:3,[0,D_coef,0])
hold on
errorbar(2,D_coef,D_std1)
text(3,D_coef,num2str(m),'fontsize',18)
ylabel('Diffusion coefficient(\mum2/s)')
saveas(gcf,[folder_name '/' bf_folder '/zbf_hist.tiff'])



end

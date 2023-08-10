function PALM2Dtracking_combine8(folder,first,last,exposure_time,trajectory_length)
close
exp_time=exposure_time; %% second
tj_len=trajectory_length-1;
%% combine all the data 
MSD_matrixb=[];
bf_folder='bfbefore';
a=regexp(folder,'\D');
folder1=folder(1:a(end));
d_sum=[];
for i=first:1:last
folder_name=[folder1 num2str(i)];  
d_sum1=load([folder_name '/' bf_folder '/zdiffusion_coefficient_individual.txt']);
d_sum=[d_sum;d_sum1];
MSD_matrixb1=load([folder_name '/' bf_folder '/' 'zMSD_matrixb.txt']);
[m,n]=size(MSD_matrixb1);
MSD_matrixb2=padarray(MSD_matrixb1,[0 tj_len-n],'post');
MSD_matrixb=[MSD_matrixb;MSD_matrixb2];
end
save([folder_name '/bf_MSD_matrixb.txt'],'-ASCII','-TABS','MSD_matrixb');
save([folder_name '/bf_d_individual.txt'],'-ASCII','-TABS','d_sum');
%% brownian
[m,n]=size(MSD_matrixb);
MSD_ave=[];
MSD_std=[];
for i=1:1:n
MSD_col= MSD_matrixb(:,i);
[j,k]=find(MSD_col>0);
if isempty(j)
    break
else
p=1:1:length(j);
MSD_ave(i)=mean(MSD_col(j(p)));
MSD_std(i)=std(MSD_col(j(p)));
end
end
MSD_aveb=MSD_ave';
MSD_stdb=MSD_std';
save([folder_name  '/bf_MSD_aveb.txt'],'-ASCII','-TABS','MSD_aveb');
save([folder_name  '/bf_MSD_stdb.txt'],'-ASCII','-TABS','MSD_stdb');
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
save([folder_name  '/bf_MSD_matrixb_bootstrap.txt'],'-ASCII','-TABS','MSD_bootstrap');

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
save([folder_name  '/bf_diffusion_coefficient.txt'],'-ASCII','-TABS','D_sum');
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
saveas(gcf,[folder_name '/bf_hist.tiff'])


end


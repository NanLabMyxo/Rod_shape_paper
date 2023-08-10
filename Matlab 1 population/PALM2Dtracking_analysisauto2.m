%% one population
clc
clear
first_frame=10;
last_frame=100;
values=0; %% 1 to view the plot, 1 not to view the plot
%%
pixel_size=0.16;%% um
exposure_time=0.1;
trajectory_length_low=4;
trajectory_length_high=12;
step_long=2;
step_short=0;
first_movie=1;
% last_movie=20;
%%
key_word='20180223mray100ms2';
folder_name='/Volumes/Elements/Guo Fu Data/PALM Nan lab/MraY/test';
d3=dir(folder_name);
for i=3:1:length(d3)
    if contains(d3(i).name,key_word)==1
       sub_folder=d3(i).name;
       folder=[folder_name '/' sub_folder]
%% select cells
PALM2Dtracking_cellselection4(folder)  
%% select outline
PALM2Dtracking_celloutline3(folder)    
    end
end

%%
for i=3:1:length(d3)
    if contains(d3(i).name,key_word)==1 
       sub_folder=d3(i).name;
       folder=[folder_name '/' sub_folder];
       d4=dir(folder);
       pp=0;
       for ii=1:1:length(d4)
           if contains(d4(ii).name,'bf')==0 && contains(d4(ii).name,'DS_Store')==0 
               pp=pp+1;
           end
       end
     last_movie=pp-2
%% find spots
PALM2Dtracking_spotfinding4(folder,first_frame,last_frame,pixel_size,values) 
%% analyze
PALM2Dtracking_sortingautomatically9(folder,pixel_size,exposure_time,trajectory_length_low,trajectory_length_high,step_long,step_short,values,first_movie,last_movie)
    end
end

%% combine data
k=1;
directory=struct([]);
for i=3:1:length(d3)
    if contains(d3(i).name,key_word)==1
    sub_folder=d3(i).name
directory(k).name=sub_folder(1:end);
k=k+1;
    end
end


%% 

trajectory_length=12;
d=directory(1).name; %% first folder
b=regexp(d(end-1:end),'\d'); %% last two digit number
n=length(d)-length(b);
% length(directory)
if length(directory)==1
    if length(b)==1
    first=d(end);
    last=d(end);
    end
    if length(b)==2
    first=d(end-1,end);
    last=d(end-1,end);
    end
sub_folder=directory(1).name;
folder=[folder_name '/' sub_folder];
PALM2Dtracking_combine8(folder,first,last,exposure_time,trajectory_length);
end
first=1;
w=1;
for i=1:1:length(directory)-1
%     i
if strncmp(directory(i).name,directory(i+1).name,n)==1
w=w+1;
if i==length(directory)-1
last=w
sub_folder=directory(i).name;
folder=[folder_name '/' sub_folder];
PALM2Dtracking_combine8(folder,first,last,exposure_time,trajectory_length);
end
else
last=w
sub_folder=directory(i).name;
folder=[folder_name,'/' sub_folder]
pause(3)
PALM2Dtracking_combine8(folder,first,last,exposure_time,trajectory_length);
w=1;
d=directory(i+1).name;
b=regexp(d(end-1:end),'\d');
n=length(d)-length(b);
end %% strncmp

end %% for

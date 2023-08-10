function cellLWratio4combine(folder,first,last)
a=regexp(folder,'\D');
folder1=folder(1:a(end));
data_sum=[];
for i=first:1:last
folder_name=[folder1 num2str(i)];
data=load([folder_name '/data.txt']);
data_sum=[data_sum;data];
end
save([folder_name '/datacombined.txt'],'-ASCII','-TABS','data_sum')
end




clear;
count=1;
num=        1213;
fig=figure;
frames(num+1)=struct('cdata',[],'colormap',[]);
for i=0:num
if i>count*num/10
message=[sprintf('%2u0',count),'% is completed'];
disp(message);
count=count+1;
end
filenum=sprintf('%04u0',i);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig\sig_',filenum);
data=load(filename);
[sigmax,index]=max(data(:,2));
plot(data(:,1),data(:,2),data(index,1),data(index,2),'pentagram');
xlim([0       0.100000000000E+00]);
ylim([     -0.583363094479E-04       0.583363094479E-04]);
drawnow;
frames(i+1)=getframe(fig);
end
video=VideoWriter('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig.mp4','MPEG-4');
video.Quality=100;
video.FrameRate=160;
open(video);
writeVideo(video,frames);
close(video);
disp('sig.mp4 is created');
clear;
count=1;
num=        1213;
fig=figure;
frames(num+1)=struct('cdata',[],'colormap',[]);
for i=0:num
if i>count*num/10
message=[sprintf('%2u0',count),'% is completed'];
disp(message);
count=count+1;
end
filenum=sprintf('%04u0',i);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\c\c_',filenum);
data=load(filename);
plot(data(:,1),data(:,2));
xlim([0       0.100000000000E+00]);
ylim([-0.01 1.01]);
drawnow;
frames(i+1)=getframe(fig);
end
video=VideoWriter('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\c.mp4','MPEG-4');
video.Quality=100;
video.FrameRate=160;
open(video);
writeVideo(video,frames);
close(video);
disp('c.mp4 is created');
clear;
fig=figure;
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\energy');
data_1=load(filename);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\energy_total');
data_2=load(filename);
plot(data_1(:,1),data_1(:,2),data_1(:,1),data_1(:,3),data_1(:,1),data_1(:,4),data_1(:,1),data_2(:,1));
xlim([0       0.636982500000E-02]);
legend('kinetic','strain','fracture','total','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\energy.jpg','-r600');
disp('energy.jpg is created');
clear;
count=1;
n_3=           0;
n_4=           0;
n_5=           0;
num_2=           0;
num_3=          25;
num_4=          50;
num_5=          75;
num_6=         100;
fig=figure;
for i=0:n_5
filenum=sprintf('%04u0',i);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\u\u_',filenum);
data=load(filename);
sig(i+1,1)=data(          51,2);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig\sig_',filenum);
data=load(filename);
sig(i+1,2)=data(num_2,2);
sig(i+1,3)=data(num_3,2);
sig(i+1,3)=(sig(i+1,3)+data(num_3+1,2))/2;
sig(i+1,4)=data(num_4,2);
sig(i+1,4)=(sig(i+1,4)+data(num_4+1,2))/2;
sig(i+1,5)=data(num_5,2);
sig(i+1,5)=(sig(i+1,5)+data(num_5+1,2))/2;
sig(i+1,6)=data(num_6,2);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\c\c_',filenum);
data=load(filename);
c(i+1,2)=data(           0,2);
c(i+1,3)=data(          13,2);
c(i+1,4)=data(          26,2);
c(i+1,5)=data(          38,2);
c(i+1,6)=data(          51,2);
end
plot(sig(1:n_3+1,1),sig(1:n_3+1,2),sig(1:n_3+1,1),sig(1:n_3+1,3),sig(1:n_3+1,1),sig(1:n_3+1,4),sig(1:n_3+1,1),sig(1:n_3+1,5),sig(1:n_3+1,1),sig(1:n_3+1,6));
xlim([0 sig(n_3+1,1)]);
ylim([0       0.116672618896E+07]);
xline(      0.530330085890E-05,'--');
yline(      0.106066017178E+07,'--');
legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig-u_00000
.jpg','-r600');
plot(sig(1:n_4+1,1),sig(1:n_4+1,2),sig(1:n_4+1,1),sig(1:n_4+1,3),sig(1:n_4+1,1),sig(1:n_4+1,4),sig(1:n_4+1,1),sig(1:n_4+1,5),sig(1:n_4+1,1),sig(1:n_4+1,6));
xlim([0 sig(n_4+1,1)]);
ylim([0       0.116672618896E+07]);
xline(      0.530330085890E-05,'--');
yline(      0.106066017178E+07,'--');
legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig-u_00000
.jpg','-r600');
plot(sig(1:n_5+1,1),sig(1:n_5+1,2),sig(1:n_5+1,1),sig(1:n_5+1,3),sig(1:n_5+1,1),sig(1:n_5+1,4),sig(1:n_5+1,1),sig(1:n_5+1,5),sig(1:n_5+1,1),sig(1:n_5+1,6));
xlim([0 sig(n_5+1,1)]);
ylim([0       0.116672618896E+07]);
xline(      0.530330085890E-05,'--');
yline(      0.106066017178E+07,'--');
legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig-u_00000
.jpg','-r600');
disp('sig-u.jpg is created');
plot(sig(1:n_3+1,2),c(1:n_3+1,2),sig(1:n_3+1,3),c(1:n_3+1,3),sig(1:n_3+1,4),c(1:n_3+1,4),sig(1:n_3+1,5),sig(1:n_3+1,5),sig(1:n_3+1,1),sig(1:n_3+1,6));
xlim([0       0.116672618896E+07]);
xline(      0.106066017178E+07,'--');
yline(0.75,'--');
legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig-u_00000
.jpg','-r600');
plot(sig(1:n_4+1,2),c(1:n_4+1,2),sig(1:n_4+1,3),c(1:n_4+1,3),sig(1:n_4+1,4),c(1:n_4+1,4),sig(1:n_4+1,5),sig(1:n_4+1,5),sig(1:n_4+1,1),sig(1:n_4+1,6));
xlim([0       0.116672618896E+07]);
xline(      0.106066017178E+07,'--');
yline(0.75,'--');
legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig-u_00000
.jpg','-r600');
plot(sig(1:n_5+1,2),c(1:n_5+1,2),sig(1:n_5+1,3),c(1:n_5+1,3),sig(1:n_5+1,4),c(1:n_5+1,4),sig(1:n_5+1,5),sig(1:n_5+1,5),sig(1:n_5+1,1),sig(1:n_5+1,6));
xlim([0       0.116672618896E+07]);
xline(      0.106066017178E+07,'--');
yline(0.75,'--');
legend('x=0','x=L/4','x=L/2','x=3L/4','x=L','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000','#FFF100'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.28.16.49\sig-u_00000
.jpg','-r600');
disp('c-sig.jpg is created');

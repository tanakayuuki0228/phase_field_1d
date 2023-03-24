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
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\sig\sig_',filenum);
data=load(filename);
[sigmax,index]=max(data(:,2));
plot(data(:,1),data(:,2),data(index,1),data(index,2),'pentagram');
xlim([0       0.100000000000E+00]);
ylim([     -0.116672618896E+07       0.116672618896E+07]);
drawnow;
frames(i+1)=getframe(fig);
end
video=VideoWriter('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\sig.mp4','MPEG-4');
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
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\c\c_',filenum);
data=load(filename);
plot(data(:,1),data(:,2));
xlim([0       0.100000000000E+00]);
ylim([-0.01 1.01]);
drawnow;
frames(i+1)=getframe(fig);
end
video=VideoWriter('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\c.mp4','MPEG-4');
video.Quality=100;
video.FrameRate=160;
open(video);
writeVideo(video,frames);
close(video);
disp('c.mp4 is created');
clear;
fig=figure;
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\energy');
data_1=load(filename);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\energy_total');
data_2=load(filename);
plot(data_1(:,1),data_1(:,2),data_1(:,1),data_1(:,3),data_1(:,1),data_1(:,4),data_1(:,1),data_2(:,1));
xlim([0       0.636982500000E-02]);
legend('kinetic','strain','fracture','total','Location','NorthEastOutside');
newcolors={'#FF4B00','#005AFF','#03AF7A','#000000'};
colororder(newcolors);
drawnow;
print(fig,'-djpeg','C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.24.21.01\energy.jpg','-r600');
disp('energy.jpg is created');

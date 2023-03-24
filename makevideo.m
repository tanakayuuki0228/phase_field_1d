clear;
count=1;
num=          10;
fig=figure;
frames(num+1)=struct('cdata',[],'colormap',[]);
for i=0:num
if i>count*num/10
message=[sprintf('%2u0',count),'% is completed'];
disp(message);
count=count+1;
end
filenum=sprintf('%04u0',i);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.22.22.33\sig\sig_',filenum);
data=load(filename);
[sigmax,index]=max(data(:,2));
plot(data(:,1),data(:,2),data(index,1),data(index,2),'pentagram');
xlim([0       0.100000000000E+00]);
ylim([     -0.116672618896E+07       0.116672618896E+07]);
drawnow;
frames(i+1)=getframe(fig);
end
video=VideoWriter('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.22.22.33\sig.mp4','MPEG-4');
video.Quality=100;
video.FrameRate=160;
open(video);
writeVideo(video,frames);
close(video);
disp('sig.mp4 is created');
clear;
count=1;
num=          10;
fig=figure;
frames(num+1)=struct('cdata',[],'colormap',[]);
for i=0:num
if i>count*num/10
message=[sprintf('%2u0',count),'% is completed'];
disp(message);
count=count+1;
end
filenum=sprintf('%04u0',i);
filename=append('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.22.22.33\c\c_',filenum);
data=load(filename);
plot(data(:,1),data(:,2));
xlim([0       0.100000000000E+00]);
ylim([-0.01 1.01]);
drawnow;
frames(i+1)=getframe(fig);
end
video=VideoWriter('C:\Users\tanaka\Documents\phase-field_1d_results\2023.03.22.22.33\c.mp4','MPEG-4');
video.Quality=100;
video.FrameRate=160;
open(video);
writeVideo(video,frames);
close(video);
disp('c.mp4 is created');

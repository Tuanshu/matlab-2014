tic
clear all

%% 8-bit binary
%tic

cd('C:\Personal\TuanShu\140924_Oral_250fps');

mkdir('divide_8bit')
N_frame=60;
size_x=1000; %140
size_y=1000; %400
num_of_division=10;

for slice_mun=1:N_frame
fin=fopen([' Oral_',sprintf('%d',slice_mun),'.']);
A=fread(fin,[size_x,size_y*num_of_division],'uint8','b');
%imagesc(A);
colormap gray

for i=1:num_of_division
    slice=A(:,(size_y*(i-1)+1):size_y*i);
%     dlmwrite([cd,'\divide\',sprintf('%d',(slice_mun-1)*10+i),'.txt'],slice);
    imwrite(slice,[cd,'\divide_8bit\',sprintf('%d',(slice_mun-1)*num_of_division+i),'.bmp']);
end

fclose('all');
% imagesc(slice);
% colormap gray

 end

%toc 

 
%% 16-bit binary
%tic

cd('C:\Personal\TuanShu\140924_Oral_250fps');

mkdir('divide')
N_frame=249;
size_x=648; %140
size_y=488; %400
num_of_division=1;

for slice_mun=1:N_frame
fin=fopen([' Oral_',sprintf('%d',slice_mun),'.']);
A=fread(fin,[size_x,size_y*num_of_division],'*uint16','b');
%imagesc(A);
colormap gray

for i=1:num_of_division
    slice=A(:,(size_y*(i-1)+1):size_y*i);
%     dlmwrite([cd,'\divide\',sprintf('%d',(slice_mun-1)*10+i),'.txt'],slice);
    imwrite(slice,[cd,'\divide\',sprintf('%d',(slice_mun-1)*num_of_division+i),'.png']);
end

fclose('all');
% imagesc(slice);
% colormap gray

 end

%toc 

%% binary to image for FFOCT-1
clear all
tic
ROI=[9 488 9 648];

cd('D:\Users\TuanShu\141020_CeYSO BCC\bin_CeYSO_125fps_25fpp_1micronsecond_log');

mkdir('divide')
N_frame=2190;

% fin=fopen([sprintf('%08d',(1000)),'.bin']);
%  A=(fread(fin,[648 488],'uint16','b'))';
%MAX=0;
%for slice_num=1:N_frame
% fin=fopen([sprintf('%08d',(slice_num)),'.bin']);
% A=(fread(fin,[648 488],'*uint32','b'))';
% B=A(ROI(1):ROI(2),ROI(3):ROI(4));
% if max(max(B))>MAX
%     MAX=max(max(B));
% end
% fclose('all');
% disp(slice_num);
%end



for slice_num=1:N_frame
 fin=fopen([sprintf('%08d',(slice_num)),'.bin']);
 A=(fread(fin,[648 488],'uint32','b'))';
 imwrite(A(ROI(1):ROI(2),ROI(3):ROI(4))/(2^32),[cd,'\divide\',sprintf('%d',slice_num),'.png'],'BitDepth',16);
 fclose('all');
end

imagesc(A(ROI(1):ROI(2),ROI(3):ROI(4)));
caxis([1E9 1.1E9]);
toc 

%% binary to image for FFOCT-1  (for noncontinuous binary file)
clear all
tic

cd('C:\Personal\Yao sheng\(©s®x)tissue carrier data\2.carrier_PZT not move1\carrier data');

mkdir('divide')
N_frame=5400;
% 
% [fin,msg]=fopen([sprintf('%08d',(67)),'.bin']);
%  A=(fread(fin,[648 488],'uint16','b'))';

 k=0;
for slice_num=1:N_frame
 fin=fopen([sprintf('%08d',(slice_num)),'.bin']);
 if fin ==-1
     k=k+1;
     fclose('all');
 else
     A=(fread(fin,[648 488],'uint16','b'))';
     imwrite(A/65535,[cd,'\divide\',sprintf('%d',slice_num-k),'.png'],'BitDepth',16);
     fclose('all');
 end
end

toc 
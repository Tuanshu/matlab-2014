tic
clear all

%% 8-bit binary
N=50;
path=sprintf('D:\\Users\\TuanShu\\140626_RBC\\carrier data\\%08d.bin',N);
%path='D:\Users\TuanShu\140626_RBC\1002.bin';

N_frame=1;
width=648;
num_of_division=10;
k=0;


for slice_num=1:N_frame
 fin=fopen(path);
 if fin ==-1
     k=k+1;
     fclose('all');
 else
     A=(fread(fin,[width inf],'uint16','b'))';
     %imwrite(A/65535,[cd,'\divide\',sprintf('%d',slice_num-k),'.png'],'BitDepth',16);
     imagesc(A);
     fclose('all');
 end
end

toc 
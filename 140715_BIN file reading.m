clear all

N=1002;
path=sprintf('D:\\Users\\TuanShu\\140626_RBC\\%d.bin',N);
fileIDD=fopen(path,'r','a');

width=648;


image=fread(fileIDD,[width,inf], 'uint16',0,'a');

imagesc(image);
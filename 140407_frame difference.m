clear all

cd('D:\Users\TuanShu\140407');

%Frame_1=dlmread('140324_fundus_2_18');   TiSa
%Frame_2=dlmread('140324_fundus_2_bg_51');   TiSa
%Frame_1=dlmread('140321_fundus_2_592');    SLD
Frame_1=dlmread('140403_Bright field_Mirror_int50_954');   
%Frame_2=dlmread('140324_fundus_2_bg_51');   
Frame_Diff=Frame_1;%-Frame_2;
Frame_Record=Frame_Diff(1:1500,1:1000);
Frame_Record=Frame_Record./max(Frame_Record(:));
subplot(1,3,1)
imagesc(Frame_1);
axis equal
xlim([0 size(Frame_1,2)]);
ylim([0 size(Frame_1,1)]);

subplot(1,3,2)

imagesc(Frame_2);
axis equal
xlim([0 size(Frame_1,2)]);
ylim([0 size(Frame_1,1)]);


subplot(1,3,3)

imagesc(Frame_Diff);
axis equal
xlim([0 size(Frame_1,2)]);
ylim([0 size(Frame_1,1)]);
caxis([0 1000]);

figure;


imagesc(Frame_Record);
axis equal
xlim([0 size(Frame_Record,2)]);
ylim([0 size(Frame_Record,1)]);
caxis([0 1]);
colormap(gray);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'XColor','white');
set(gca,'YColor','white');
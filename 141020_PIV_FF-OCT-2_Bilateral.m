clear all

%% Note: ¶Âªº¤~¬Ocell

interrogation_size=20;

Total_Size_X=488;
Total_Size_Y=648;

Pixel_Spacing=0.45;     %micron

Frame_Rate=250;%*526/2000;         %fps
Normal_Velocity=Pixel_Spacing*Frame_Rate;  %micron/second
Path='D:\Users\TuanShu\141020_CeYSO BCC\bin_CeYSO_125fps_25fpp_1micronsecond_log\';
cd(Path);
File_list=dir;
%MOVIE=aviread('red blood cell in tissue.avi');

%Test=MOVIE(205).cdata(:,:,1)-MOVIE(205).cdata(:,:,3);

%%
ROI=[9 488 9 648];

N_frame=1759;



for slice_num=1:N_frame
 fin=fopen([sprintf('%08d',(slice_num)),'.bin']);
 A=(fread(fin,[648 488],'uint32','b'))';
 Image(:,:,slice_num)=A(ROI(1):ROI(2),ROI(3):ROI(4));
 fclose('all');
 disp(slice_num);
end



%%
%LP=3;
%FFT_Image=fft(Image,[],3);

%FFT_Image(:,:,1:LP)=0;
%FFT_Image(:,:,(size(FFT_Image,3)/2+1):size(FFT_Image,3))=0;

%Image_New=abs(ifft(FFT_Image,[],3));

%%
N=928;
%Image_New(Image_New<Threshold)=0;

Image_Temp=Image(:,:,N);
       
imagesc(log10(Image_Temp));
colormap(gray);
%caxis([0 0.1]);
%%
cd('D:\Users\TuanShu\141020_CeYSO BCC\Bilateral Filtering');

W=10;
SIGMA=[5 0.03];
N=928;
%Image_New(Image_New<Threshold)=0;

        Image_Temp=Image(:,:,N)/max(max(Image(:,:,N)));
        Image_Bi_2D=bfilter2(Image_Temp,W,SIGMA);
        
subplot(2,1,1)
imagesc(Image_Temp);

subplot(2,1,2)
imagesc(log10(Image_Bi_2D));
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
%% All Bilateral
cd(Path);
Start_N=200;
End_N=750;
Image_Bi=Image(:,:,Start_N:End_N);
mkdir('\Bilateral\'); 

for p=Start_N:End_N
    Image_Bi(:,:,p)=bfilter2(Image(:,:,p),W,SIGMA);
    imwrite(Image_Bi(:,:,p),[cd,'\Bilateral\',sprintf('%d',p),'.png'],'BitDepth',16);

    disp(p);
end
%%
N=719;

subplot(2,1,1)
imagesc(Image(:,:,N));

subplot(2,1,2)
imagesc(Image_Bi(:,:,N));
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');

%%
%% backup
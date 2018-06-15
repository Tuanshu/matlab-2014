clear all

cd('D:\Users\TuanShu\131111\');

Calibration_Array=importdata('Calibration_Array.txt');
fx=importdata('fx_new.txt');

%Array_Number=119:273;
Array_Number=1:100;

LP=200;
%%
Image(1:2048,1:1088)=0;
for p=1:length(Array_Number)
    Data1_2D_Temp=importdata(sprintf('131111_paper_inter_%i',Array_Number(p)));
    %Data1_2D_Temp=importdata(sprintf('131111_onion_inter_%i',Array_Number(p)));
    Data1_2D_New=interp1(Calibration_Array,Data1_2D_Temp,fx);
    Data1_2D_New(isnan(Data1_2D_New))=0;

    %imagesc(Data1_2D_New);
    %axis equal
    %colormap(gray);
    FFT = abs(fft(Data1_2D_New,[],1));
    FFT=FFT(1:round(size(FFT,1)/2),:);
    FFT(1:LP,:)=0;
    Image=Image+FFT;
    disp(p);
end

Image=Image-min(min(Image));
Image_dB=20*log10(Image);
Image_dB=Image_dB-max(max(Image_dB));
imagesc(Image_dB,'xdata',4.65*[1:size(Image,2)],'ydata',0.2*[1:size(Image,1)]);
%axis equal
colormap(gray);

Image_single=FFT;
Image_single_dB=20*log10(Image_single);
Image_single_dB=Image_single_dB-max(max(Image_single_dB));
imagesc(Image_single_dB,'xdata',4.65*[1:size(Image,2)],'ydata',0.2*[1:size(Image,1)]);
%axis equal

plot(real(FFT_1(:,500)))
c=3E8;
Time_total=1/(max(fx)/(length(fx)-1));
Time=[0:Time_total/(length(fx)-1):Time_total]/2;%/2�O�]���@�Ӥ@�^
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;
d_Position_micron=Position_micron(2)-Position_micron(1);
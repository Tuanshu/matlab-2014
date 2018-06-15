clear all

cd('D:\Users\TuanShu\140114\');
C=3E8;

LP=1000;


SP=2;
%% 注意 原本的架構, pixel number小> frequency大, 但因為cuda的程式似乎不接受frequency由大到小的calibration array, 所以旋轉camera後重新量

M=620;  %which line for calibration
N=1000; %which pixel is center

Frequency_max=8E14;

N_f=4096;

fx=[0:Frequency_max/(N_f-1):Frequency_max]';
dlmwrite('fx_new.txt',fx,'delimiter','\t','newline','pc');

%%fx=importdata('fx_new.txt');
%4  43.7141
%5  93.7087
%6  135.5168


Thickness_1=43.7141*1E-6;
Thickness_2=93.7087*1E-6;

Data1_2D=importdata('140114_cali_4_44.txt');

Data1_2D=rot90(Data1_2D);

Data2_2D=importdata('140114_cali_5_57.txt');

Data2_2D=rot90(Data2_2D);


imagesc((Data2_2D));
colormap(gray);
xlabel('X Pixel Number');
ylabel('Y Pixel Number');
%%
Data1=Data1_2D(M,:);
Data2=Data2_2D(M,:);

X=1:length(Data1);
%%
%Thickness_1=-Thickness_1;
%Thickness_2=-Thickness_2;

plot(X,Data1,X,Data2);
%%

FFT1=fft(Data1);


FFT1(1:SP)=0;
FFT1(LP:end)=0;
Data1_H=ifft(FFT1);
PHI1=unwrap(angle(Data1_H));

FFT2=fft(Data2);
FFT2(1:SP)=0;
FFT2(LP:end)=0;
Data2_H=ifft(FFT2);
PHI2=unwrap(angle(Data2_H));

Delta_PHI=PHI2-PHI1;


plot(X,FFT1,X,FFT2);

plot(X,PHI1,X,PHI2);
xlim([0 2040]);
xlabel('X Pixel Number');
ylabel('Phase (radian)');
%%
Freq=Delta_PHI*C/(Thickness_2-Thickness_1)/4/pi;

plot(X,Delta_PHI);
plot(Freq,Delta_PHI);

Freq=Freq-Freq(N)+C/(780E-9);

plot(X,Freq);
%plot(Freq);
xlim([0 2040]);
xlabel('X Pixel Number');
ylabel('Optical Frequency (Hz)');

Calibration_array=Freq';
Calibration_array_New=smooth(Calibration_array,100);
index_array=[1:length(Calibration_array_New)]';
index_start=300;
index_end=1700;
XX=fit(index_array(index_start:index_end),Calibration_array_New(index_start:index_end),'poly2');
Calibration_array_Fit=XX.p1*index_array.^2 + XX.p2*index_array + XX.p3;

plot(index_array,Calibration_array_New,index_array,Calibration_array_Fit);
xlabel('X Pixel Number');
ylabel('Optical Frequency (Hz)');
xlim([0 2040]);

dlmwrite('Calibration_array.txt',Calibration_array_Fit,'delimiter','\t','newline','pc');

%dlmwrite('Calibration_array.txt',Calibration_array_New,'delimiter','\t','newline','pc');

plot(smooth(Calibration_array,100),Data2);
xlabel('Optical Frequency (Hz)');
ylabel('(LSB)');
%%
for p=1:size(Data1_2D,1)
    Data1_2D_New(p,:)=interp1(Calibration_array,Data1_2D(p,:),fx);
    Data2_2D_New(p,:)=interp1(Calibration_array,Data2_2D(p,:),fx);
    disp(p);
end
Data1_2D_New(isnan(Data1_2D_New))=0;
Data2_2D_New(isnan(Data2_2D_New))=0;

image(Data2_2D_New);
colormap(gray);
FFT_1 = abs(fft(Data1_2D_New,[],2));
FFT_2 = abs(fft(Data2_2D_New,[],2));

image(FFT_2);
colormap(gray);

plot(real(FFT_2(500,:)))
c=3E8;
Time_total=1/(max(fx)/(length(fx)-1));
Time=[0:Time_total/(length(fx)-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;
d_Position_micron=Position_micron(2)-Position_micron(1);
clear all

Sampling_rate=200000;
Background_Smooth_Window=500;


cd('D:\Users\TuanShu\140604\exp5\');

Sets=50;

for p=1:Sets
    Data_1(:,p)=dlmread(sprintf('Set_%02d_Profile_1.txt',p)); 
    Data_2(:,p)=dlmread(sprintf('Set_%02d_Profile_2.txt',p)); 
    Data_3(:,p)=dlmread(sprintf('Set_%02d_Profile_3.txt',p)); 
    Data_1(:,p)=Data_1(:,p)-Data_1(1,p); 
    Data_2(:,p)=Data_2(:,p)-Data_2(1,p); 
    Data_3(:,p)=Data_3(:,p)-Data_3(1,p); 
end

Time=[1:size(Data_1,1)]/Sampling_rate*1000;     %ms

plot(Time,Data_1);

%%
N=9;
plot(Time,Data_1(:,N));

%%
LowerBand=25000;    %(Hz)
UpperBand=35000;    %(Hz)


Frequency=[1:size(Data_1,1)]*Sampling_rate/size(Data_1,1);

FFT_Data_1=fft(Data_1,[],1);
FFT_Data_2=fft(Data_2,[],1);
FFT_Data_3=fft(Data_3,[],1);

plot(Frequency,FFT_Data_1);

FFT_Data_1_Filtered=FFT_Data_1;
FFT_Data_2_Filtered=FFT_Data_2;
FFT_Data_3_Filtered=FFT_Data_3;
FFT_Data_1_Filtered(Frequency<LowerBand,:)=0;
FFT_Data_1_Filtered(Frequency>UpperBand,:)=0;
FFT_Data_2_Filtered(Frequency<LowerBand,:)=0;
FFT_Data_2_Filtered(Frequency>UpperBand,:)=0;
FFT_Data_3_Filtered(Frequency<LowerBand,:)=0;
FFT_Data_3_Filtered(Frequency>UpperBand,:)=0;

Data_1_Filtered=ifft(FFT_Data_1_Filtered,[],1);
Data_2_Filtered=ifft(FFT_Data_2_Filtered,[],1);
Data_3_Filtered=ifft(FFT_Data_3_Filtered,[],1);

%% Smoothing
Smooth_Window=50;
for p=1:Sets
    Data_1_Filtered_Envelope(:,p)=smooth(abs(Data_1_Filtered(:,p)),Smooth_Window,'sgolay');
    Data_2_Filtered_Envelope(:,p)=smooth(abs(Data_2_Filtered(:,p)),Smooth_Window,'sgolay');
    Data_3_Filtered_Envelope(:,p)=smooth(abs(Data_3_Filtered(:,p)),Smooth_Window,'sgolay');

end
%%
N=9;
subplot(3,1,1)
plot(Time,Data_1_Filtered_Envelope(:,N),Time,Data_1_Filtered(:,N));
ylim([0 15E-3]);

subplot(3,1,2)
plot(Time,Data_2_Filtered_Envelope(:,N),Time,Data_2_Filtered(:,N));
ylim([0 15E-3]);

subplot(3,1,3)
plot(Time,Data_3_Filtered_Envelope(:,N),Time,Data_3_Filtered(:,N));
ylim([0 15E-3]);

%%
N=9;
subplot(2,1,1)
plot(Time,Data_1(:,N));

subplot(2,1,2)

plot(Time,Data_1_Filtered_Envelope(:,N),Time,Data_1_Filtered(:,N));
ylim([0 10E-3]);


cd('D:\Users\TuanShu\');
dlmwrite('Time.txt',Time','delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite('Data_1.txt',Data_1(:,N),'delimiter','\t','newline','pc','precision', '%.12f');
dlmwrite('Data_1_Filtered_Envelope.txt',Data_1_Filtered_Envelope(:,N),'delimiter','\t','newline','pc','precision', '%.12f');

%%
N=17;
plot(Time,Data_2_Filtered(:,:));

%% Peak finding

MINPEAKHEIGHT=0.5E-3;
Range_Position=10;
[value Range]=find(Time>Range_Position,1,'first');
Starting_Position_1=16;
Starting_Position_2=15;
Starting_Position_3=12;

Starting_Index_1=find(Time>Starting_Position_1,1,'first');
Starting_Index_2=find(Time>Starting_Position_2,1,'first');
Starting_Index_3=find(Time>Starting_Position_3,1,'first');


for p=2:Sets
    if p==2
        if isempty(findpeaks(Data_1_Filtered_Envelope(Starting_Index_1:(Starting_Index_1+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT))~=1
            [value Data_1_Inward_Index(p)]=findpeaks(Data_1_Filtered_Envelope(Starting_Index_1:(Starting_Index_1+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT);
            Data_1_Inward_Index(p)=Data_1_Inward_Index(p)+Starting_Index_1-1;
        else
            Data_1_Inward_Index(p)=1;
        end
        
        if isempty(findpeaks(Data_2_Filtered_Envelope(Starting_Index_2:(Starting_Index_2+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT))~=1
            [value Data_2_Inward_Index(p)]=findpeaks(Data_2_Filtered_Envelope(Starting_Index_2:(Starting_Index_2+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT);
            Data_2_Inward_Index(p)=Data_2_Inward_Index(p)+Starting_Index_2-1;
        else
            Data_2_Inward_Index(p)=1;
        end
        
        if isempty(findpeaks(Data_3_Filtered_Envelope(Starting_Index_3:(Starting_Index_3+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT))~=1
            [value Data_3_Inward_Index(p)]=findpeaks(Data_3_Filtered_Envelope(Starting_Index_3:(Starting_Index_3+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT);
            Data_3_Inward_Index(p)=Data_3_Inward_Index(p)+Starting_Index_3-1;

        else
            Data_3_Inward_Index(p)=1;
        end
    else
        if (isempty(findpeaks(Data_1_Filtered_Envelope(Data_1_Inward_Index(p-1):(Data_1_Inward_Index(p-1)+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT))~=1)&&(Data_1_Inward_Index(p-1)+Range)<size(Data_1_Inward_Index,1)
            [value Data_1_Inward_Index(p)]=findpeaks(Data_1_Filtered_Envelope(Data_1_Inward_Index(p-1):(Data_1_Inward_Index(p-1)+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT);
            Data_1_Inward_Index(p)=Data_1_Inward_Index(p)+Data_1_Inward_Index(p-1)-1;
        else
            Data_1_Inward_Index(p)=Data_1_Inward_Index(p-1);
        end
        
        if (isempty(findpeaks(Data_2_Filtered_Envelope(Data_2_Inward_Index(p-1):(Data_2_Inward_Index(p-1)+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT))~=1)&&(Data_2_Inward_Index(p-1)+Range)<size(Data_2_Inward_Index,1)
            [value Data_2_Inward_Index(p)]=findpeaks(Data_2_Filtered_Envelope(Data_2_Inward_Index(p-1):(Data_2_Inward_Index(p-1)+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT);
            Data_2_Inward_Index(p)=Data_2_Inward_Index(p)+Data_2_Inward_Index(p-1)-1;
        else
            Data_2_Inward_Index(p)=Data_2_Inward_Index(p-1);
        end
        
        if (isempty(findpeaks(Data_3_Filtered_Envelope(Data_3_Inward_Index(p-1):(Data_3_Inward_Index(p-1)+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT))~=1)&&(Data_3_Inward_Index(p-1)+Range)<size(Data_3_Inward_Index,1)
            [value Data_3_Inward_Index(p)]=findpeaks(Data_3_Filtered_Envelope(Data_3_Inward_Index(p-1):(Data_3_Inward_Index(p-1)+Range),p),'NPEAKS',1,'MINPEAKHEIGHT',MINPEAKHEIGHT);
            Data_3_Inward_Index(p)=Data_3_Inward_Index(p)+Data_3_Inward_Index(p-1)-1;

        else
            Data_3_Inward_Index(p)=Data_3_Inward_Index(p-1);
        end
    end
end

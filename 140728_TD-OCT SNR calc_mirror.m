clear all

cd('D:\Users\TuanShu\140728_FF-OCT characteristic\');

Data=dlmread('140728_test_somehow optimized_1-5V_600mA_60dB.txt');        

%% Raw Data Reading

sampling_rate=128000/4000;  %Hz

scanning_speed=1;   %micron/sec

Pre_ave=1;

for p=1:length(Data)/Pre_ave
    Data_Pre(p)=sum(Data((Pre_ave*(p-1)+1):(Pre_ave*p)));
end
Data_Pre=Data_Pre';
dPosition=scanning_speed/sampling_rate*Pre_ave; %micron


X=(dPosition:dPosition:dPosition*length(Data_Pre))';


plot(X,Data_Pre);
xlabel('Position (micron)');
ylabel('Intensity (a.u.)');


%% Background Substration
Background_start=Data_Pre(1);
Background_end=Data_Pre(end);

Background=interp1([X(1) X(end)]',[Background_start Background_end]',X);

plot(X,Data_Pre,X,Background);

Data_sub=Data_Pre-Background;

plot(X,Data_sub);

%% Filtering
LB=310;
HB=500;

Data_FFT=fft(Data_sub);

plot(real(Data_FFT));

Data_FFT_Filtered=Data_FFT;

Data_FFT_Filtered(1:LB)=0;
Data_FFT_Filtered(HB:end)=0;

Data_New=ifft(Data_FFT_Filtered);

plot(X,abs(Data_New));

%% Averaging

Averging_Factor=1;

Data_New_Averaged=smooth(abs(Data_New),Averging_Factor);

% SNR calc
Noise_Floor=abs(Data_New_Averaged((length(Data_New_Averaged)-999):length(Data_New_Averaged)));

plot(Noise_Floor);
Error=std(Noise_Floor);

SNR=log10(max(abs(Data_New_Averaged))/Error)*20;


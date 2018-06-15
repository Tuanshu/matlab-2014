clear all

cd('D:\Users\TuanShu\140728_FF-OCT characteristic\');
for N=1:4
Data=dlmread(sprintf('140728_test_somehow optimized_1-5V_600mA_60dB_%d.txt',N))';        
Data=Data(1:4420);
%% Raw Data Reading

sampling_rate=128000/4000;  %Hz

scanning_speed=1;   %micron/sec

Pre_ave=1;

for p=1:length(Data)/Pre_ave
    Data_Pre(p)=sum(Data((Pre_ave*(p-1)+1):(Pre_ave*p)));
end

dPosition=scanning_speed/sampling_rate*Pre_ave; %micron


X=(dPosition:dPosition:dPosition*length(Data_Pre));


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
LB=10;
HB=1000;

Data_FFT=fft(Data_sub);

plot(real(Data_FFT));

Data_FFT_Filtered=Data_FFT;

Data_FFT_Filtered(1:LB)=0;
Data_FFT_Filtered(HB:end)=0;

Data_New(:,N)=ifft(Data_FFT_Filtered);
end



plot(X,abs(Data_New));

%% Circshift
    [Max_1 Max_1_index]=max(Data_New(:,1));
Max_1_index=1349;
for N=1:4
    [Max_N Max_N_index]=max(Data_New(:,N));
    Data_New_Normalized(:,N)=circshift(Data_New(:,N),Max_1_index-Max_N_index)/Max_N;
end

Data_New_Normalized_Mean=mean(Data_New_Normalized,2);
plot(X,abs(Data_New_Normalized));
plot(X,log10(abs(Data_New_Normalized))*20);
xlabel('Position (micron)');
ylabel('Intensity (dB)');
xlim([0 80]);
plot(X,log10(abs(Data_New_Normalized_Mean))*20);
xlabel('Position (micron)');
ylabel('Intensity (dB)');
xlim([0 80]);
%% Averaging

Averging_Factor=1;

Data_New_Averaged=smooth(abs(Data_New),Averging_Factor);

% SNR calc
Noise_Floor=abs(Data_New_Averaged((length(Data_New_Averaged)-999):length(Data_New_Averaged)));

plot(Noise_Floor);
Error=std(Noise_Floor);

SNR=log10(max(abs(Data_New_Averaged))/Error)*20;
dlmwrite('Data_New_Normalized_Mean_log.txt',log10(abs(Data_New_Normalized_Mean))*20,'delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite('X.txt',X','delimiter','\t','newline','pc','precision', '%.6f');


clear all
cd('D:\Users\TuanShu\140307\Best\');
Data=dlmread('Waveform_SR1000Hz_V1micronsec_Order2_withbackwalking.txt');
Time_ori=Data(:,1);
Voltage_ori=Data(:,2);
plot(Time_ori,Voltage_ori);

%% Sampling Rate Reduection
%Ratio_resample=10;

%Time_resample=resample(Time_ori,1,Ratio_resample);
%Voltage_resample=resample(Voltage_ori,1,Ratio_resample);


%plot(Time_ori,Voltage_ori,Time_resample,Voltage_resample);
%%
%End_Point=96;
%Time=Time_ori(Time_ori<End_Point);
%Voltage=Voltage_ori(Time_ori<End_Point);
%plot(Time,Voltage);

%plot(Voltage);
%%
Smooth_Window=100;
Voltage_Smooth=smooth(Voltage_ori,Smooth_Window,'loess');
Residual=Voltage_Smooth-Voltage_ori;
plot(Time_ori,Voltage_Smooth);

plot(Time_ori,Residual);

Output=([Time_ori Voltage_Smooth]);

dlmwrite('Waveform_SR1000Hz_V1micronsec_Order2_withbackwalking_Smooth.txt',Output,'delimiter','\t','newline','pc','precision', '%.6f');



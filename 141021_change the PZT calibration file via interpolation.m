clear all
cd('D:\Users\TuanShu\140725_PZT calibration');

%%
Data=dlmread('Waveform_SR1000Hz_V1micronsec_Order3_withbackwalking_NEW_smoothed.txt');        
plot(Data(:,1),Data(:,2));

Original_Velocity=1;

Time=Data(:,1);
Voltage=Data(:,2);

%%  filtering the signal
Ratio=0.2;      %i.e. 如果為0.5, 那就是減速一半

D_Time=Data(2,1)-Data(1,1);
Time_Max_Old=max(Time);
Time_Max_New=round(Time_Max_Old/Ratio);

Time_New_Waveform=D_Time:D_Time:Time_Max_New;

Voltage_New_Waveform=interp1(Time,Voltage,Time_New_Waveform*Ratio,'linear','extrap');

Voltage_New_Waveform(Voltage_New_Waveform<-1.999)=-1.999;

plot(Time_New_Waveform,Voltage_New_Waveform,'linewidth',2);
xlabel('Time (second)','fontsize',12);
ylabel('Voltage (V)','fontsize',12);

Output=[Time_New_Waveform' Voltage_New_Waveform'];

dlmwrite(sprintf('Waveform_SR%dHz_V%0.1dmicronsec_INTERPOLED_withbackwalking.txt',1/D_Time,Original_Velocity*Ratio),Output,'delimiter','\t','newline','pc','precision', '%.6f');

clear all
cd('D:\Users\TuanShu\140307\Best\Cut\');
Data=dlmread('Waveform_SR1000Hz_V8micronsec_Order3_withbackwalking_Smooth.txt');
Time_ori=Data(:,1);
Voltage_ori=Data(:,2);
plot(Time_ori,Voltage_ori);

%%
Cut_point_Rough=150;
Range=100;

Time_cut_rough=Time_ori(Time_ori<Cut_point_Rough);
Voltage_cut_rough=Voltage_ori(Time_ori<Cut_point_Rough);
plot(Time_cut_rough,Voltage_cut_rough);

Voltage_temp=Voltage_cut_rough;

Voltage_temp(Time_cut_rough<(Cut_point_Rough-Range))=999;

[minvalue minindex]=min(Voltage_temp);

Time_cut=Time_cut_rough(1:minindex);
Voltage_cut=Voltage_cut_rough(1:minindex);

plot(Time_cut,Voltage_cut);


Output=([Time_cut Voltage_cut]);

dlmwrite('Waveform_SR1000Hz_V8micronsec_Order3_withbackwalking_Smooth_Cut.txt',Output,'delimiter','\t','newline','pc','precision', '%.6f');



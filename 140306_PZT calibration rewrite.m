clear all
cd('D:\Users\TuanShu\140307\');
%%
Sampling_rate=1000;         %points/sec


%%
Wavelength_Laser=0.6328;     %micron

%%

SPF=5000;
LPF=9000;

Start_Pixel=400;
End_Pixel=278000;

Start_Voltage=-2;
Range=12;

End_Voltage=Start_Voltage+Range;

Smooth_Window=10;

Voltage_cutoff=0; %how many voltage change near the turning points is neglected.

%%
Data=dlmread('140307_3rd_measurement (8micronsec)_SR1000.txt');        
plot(Data(:,2));
Voltage_read=Data(Start_Pixel:End_Pixel,2);
Signal_read=Data(Start_Pixel:End_Pixel,1);
plot(Voltage_read);
%%  filtering the signal
FFT_Signal_read=fft(Signal_read);
plot(real(FFT_Signal_read));
FFT_Signal_read(1:SPF)=0;
FFT_Signal_read(LPF:end)=0;
Signal_read_new=ifft(FFT_Signal_read);

%% to delete elements from array
Voltage=smooth(Voltage_read,Smooth_Window);
Signal=smooth(Signal_read_new,Smooth_Window);
Voltage_for_turningpoint=smooth(Voltage_read,Smooth_Window*20); %*10 for 4 micron per sec, *20 for 2 micron per second

%Signal(Voltage>(End_Voltage-Voltage_cutoff))=[];    %�`�N����
%Voltage(Voltage>(End_Voltage-Voltage_cutoff))=[];

%% phase wunwrapping
Phase_waveform_original=unwrap(angle(Signal));

%% �h��turning points
Voltage_invert=max(Voltage_for_turningpoint)-Voltage_for_turningpoint;
[minvalue minindex]=findpeaks(Voltage_invert,'MINPEAKDISTANCE',10000,'MINPEAKHEIGHT',End_Voltage-1);
[maxvalue maxindex]=findpeaks(Voltage_for_turningpoint,'MINPEAKDISTANCE',10000,'MINPEAKHEIGHT',End_Voltage-1);
index_all_turningpoints=sort([minindex; maxindex]);

%% to turn the direction of phase
Phase_waveform=Phase_waveform_original;
for p=1:length(index_all_turningpoints)
    Phase_waveform((index_all_turningpoints(p)+1):end)=2*Phase_waveform((index_all_turningpoints(p)))-(Phase_waveform((index_all_turningpoints(p)+1):end));
end

%% Position

Position_waveform=-1*Wavelength_Laser/2*Phase_waveform/2/pi;
Time_waveform=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_waveform))';

plot(Position_waveform);
plot(Voltage,Position_waveform);
plot(Time_waveform,Position_waveform);


%% Calibration
Order=3;
Velocity=8;                 %micron/sec
Position_Sampling_Resolution=Velocity/Sampling_rate;

%% Forward
Nth_Start_Turning_Point=1;
Nth_End_Turning_Point=2;

Position_Patial=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Voltage_Patial=Voltage(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Position_Patial=Position_Patial(Voltage_Patial<10);
Voltage_Patial=Voltage_Patial(Voltage_Patial<10);


plot(Voltage_Patial,Position_Patial);
Position_UniformBase=(min(Position_Patial):Position_Sampling_Resolution:max(Position_Patial))';
Time_UniformBase=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_UniformBase))';
Voltage_UniformBase=interp1(Position_Patial,Voltage_Patial,Position_UniformBase);

%% Backward

Nth_Start_Turning_Point_2=2;
Nth_End_Turning_Point_2=3;

Position_Patial_2=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point_2):(index_all_turningpoints(Nth_End_Turning_Point_2)-1));
Voltage_Patial_2=Voltage(index_all_turningpoints(Nth_Start_Turning_Point_2):(index_all_turningpoints(Nth_End_Turning_Point_2)-1));
Position_Patial_2=Position_Patial_2(Voltage_Patial_2<10);
Voltage_Patial_2=Voltage_Patial_2(Voltage_Patial_2<10);

plot(Voltage_Patial,Position_Patial);
Position_UniformBase_2=(max(Position_Patial_2):(-1*Position_Sampling_Resolution):min(Position_Patial_2))';
Time_UniformBase_2=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_UniformBase_2))';
Voltage_UniformBase_2=interp1(Position_Patial_2,Voltage_Patial_2,Position_UniformBase_2);



%% to generate the 3 loops + backwalking
Number_of_Cycle=3;
Pre_Voltage=0;
%%
Voltage_half_cycle_forward=Voltage_UniformBase;
Voltage_half_cycle_backward=Voltage_UniformBase_2;
Voltage_one_cycle=[Voltage_half_cycle_forward; Voltage_half_cycle_backward];
Voltage=repmat(Voltage_one_cycle,[Number_of_Cycle 1]);

%% backwalking
Voltage_pre_cycle_backward=Voltage_half_cycle_backward(Voltage_half_cycle_backward<Pre_Voltage);

%%
Voltage_with_pre=[Voltage_pre_cycle_backward; Voltage];
Time_with_pre=(1/Sampling_rate)*(1:length(Voltage_with_pre))';

plot(Time_with_pre,Voltage_with_pre);



Output=[Time_with_pre Voltage_with_pre];

dlmwrite(sprintf('Waveform_SR%dHz_V%dmicronsec_Order%d_withbackwalking.txt',Sampling_rate,Velocity,Order),Output,'delimiter','\t','newline','pc','precision', '%.6f');


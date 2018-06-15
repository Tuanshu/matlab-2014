clear all

%% Setting
Motor_Speed=0.5;%0.0025;   %mm/sec
Integration_Time=5;   %40%total integration time,ms

%Pixel_Average_Axial=1;
%Pixel_Average_Lateral=10;

Lateral_Spacing=Motor_Speed*Integration_Time;   %micron

Max_Wavelength=1100;             %nm
Min_Wavelength=500;             %nm
N_f=8192;
N_t=8192*8;
ROI_ratio=1/8;                  %only consider the first ROI_ratio data in TD

DC_Cutoff=5;                 %micron

array=1;

cd('D:\Users\TuanShu\140114\');
%% Data Loading

Signal_Bscan_Envelope(1:round(N_t*ROI_ratio),1:length(array))=0;


%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('140114_cali_OO_6.txt');


plot(Data(:,1),Data(:,2));


Wavelength=Data(:,1);           %nm

%Data=importdata('inter.txt')-importdata('sam.txt')-importdata('ref.txt')+importdata('bs.txt');
%cd('D:\120222\');
%Data_R=importdata('R1 Ref.txt');
%cd('D:\120222\R1 Sam\');
%Data_S=importdata(sprintf('D%i.txt',array(jj)));


C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';



Spectrum=Data(:,2);
Spectrum=Spectrum-Spectrum(1);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
%plot(Frequency_New,Spectrum_New);

%% To time domain

Spectrum_New((N_f+1):N_t)=0;
Signal=fft(Spectrum_New);

Spectrum_New=Spectrum_New(1:N_f);
%Signal=downsample(conv((Signal),(ones(Pixel_Average_Axial,1))/Pixel_Average_Axial,'same'),Pixel_Average_Axial);
    Signal_DC=Signal;
    Signal_DC(1000:(length(Signal_DC)-1000))=0;

%Signal=Signal-Signal_DC;
Signal_Bscan_Envelope(:,1)=abs(Signal(1:size(Signal_Bscan_Envelope,1)));
Signal_Bscan_Carrier(:,1)=real(Signal(1:size(Signal_Bscan_Envelope,1)));

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position_micron=C*Time(1:round(length(Time)*ROI_ratio))*1E6;
Position_micron=Position_micron(1:size(Signal_Bscan_Envelope,1));

Signal_Bscan_Envelope=Signal_Bscan_Envelope/max(max(Signal_Bscan_Envelope));
Lateral_Position=[Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing]';


%%
plot(Position_micron,Signal_Bscan_Carrier(:,1));
xlabel('OPD (micron)');
ylabel('Interference Signal');
%%
%axis equal
imagesc(Signal_Bscan_Envelope,'xdata',Lateral_Spacing:Lateral_Spacing:size(Signal_Bscan_Envelope,2)*Lateral_Spacing,'ydata',Position_micron);
%colormap(gray);
%axis equal
xlabel('Lateral Position (micron)');
ylabel('Optical Path Difference (micron)');



[maxvalue maxindex]=max(Signal_Bscan_Envelope(:,1));

plot(Lateral_Position,Signal_Bscan_Envelope(maxindex,:));

%% Circshift
%imagesc((Signal_Bscan_Envelope_Shifted(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first'),:)),'xdata',Lateral_Position,'ydata',Position_micron(find(Position_micron>Range_3,1,'first'):find(Position_micron>Range_4,1,'first')));
imagesc((Signal_Bscan_Envelope),'xdata',Lateral_Position,'ydata',Position_micron);
colormap(gray);
%caxis([-20 0]);
%axis equal
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'XColor','white');
set(gca,'YColor','white');
Signal_Bscan_Envelope_Cut=Signal_Bscan_Envelope(Position_micron>DC_Cutoff,:);
Position_micron_Cut=Position_micron(Position_micron>DC_Cutoff);
[max_value max_index]=max(Signal_Bscan_Envelope_Cut);

plot(Position_micron_Cut,Signal_Bscan_Envelope_Cut(:,1));
xlabel('Optical Path Difference (micron)');
ylabel('Normalized Interference Intensity (a.u.)');
[max_value max_index]=max(Signal_Bscan_Envelope_Cut);
Position_micron_Cut(max_index)
%xlabel('Lateral Position (micron)');
%ylabel('Optical Path Difference (micron)');
%dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

%dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');
%plot(Position_micron,Signal_Bscan_Envelope);
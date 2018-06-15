clear all


Max_Wavelength=700;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*128;
ROI_ratio=1/128;                  %only consider the first ROI_ratio data in TD

DC_Cutoff=5;                 %micron

array=1500:1500;


%% Wavelength
cd('D:\Users\TuanShu\140808_CeYSO spectra\Used');
Name=[5];
for N=1:1
Data=importdata(sprintf('140805_CeYSO_smooth 3_with lens_%dmW.txt',Name(N)));

%Data=importdata('140811_test_CeYSO_spectrum_in det_ref.txt');
%Data_BS=importdata('140811_test_CeYSO_spectrum_in det_ref_bs.txt');
Wavelength=Data(:,1);          %nm

C=3E8;

Frequency=C./(Wavelength*1E-9);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position_micron=C*Time(1:round(length(Time)*ROI_ratio))*1E6;
Spectrum=Data(:,2);%-Data_BS(:,2);
Spectrum=Spectrum-min(Spectrum);
Spectrum(1:1000)=0;

plot(Wavelength,Spectrum);
%Spectrum=Spectrum-Spectrum(1);%-Data_R(:,2)-Data_S(:,2);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
Spectrum_New(Frequency_New>Max_Frequency)=0;

Wavelength_Micron=C./Frequency_New*1E6;
Spectrum_Record(:,N)=Spectrum_New;

end

plot(Wavelength_Micron(Wavelength_Micron<0.8),Spectrum_Record(Wavelength_Micron<0.8,:));
xlabel('Wavelength (micron)');
    %% To time domain

Spectrum_New((N_f+1):N_t)=0;
Signal=fftshift(fft(Spectrum_New));

Spectrum_New=Spectrum_New(1:N_f);
    %Signal=downsample(conv((Signal),(ones(Pixel_Average_Axial,1))/Pixel_Average_Axial,'same'),Pixel_Average_Axial);
    %Signal=Signal-Signal_DC;
Signal_Patial=real(Signal((1:round(N_t*ROI_ratio))+round(length(Signal)/2)-round(N_t*ROI_ratio/2)));


Signal_Patial_Envelope=abs(Signal((1:round(N_t*ROI_ratio))+round(length(Signal)/2)-round(N_t*ROI_ratio/2)));

Signal_Patial_Envelope=Signal_Patial_Envelope/max(Signal_Patial_Envelope);

plot(Position_micron-5,Signal_Patial/max(Signal_Patial),Position_micron-5,Signal_Patial_Envelope);
xlabel('Position (micron)');
ylabel('Spectral Intensity');
xlim([-5 5]);

FWHM=Position_micron(find(Signal_Patial_Envelope>0.5,1,'last'))-Position_micron(find(Signal_Patial_Envelope>0.5,1,'first'));




    clear;
    
    fid = fopen('TV_CCT.txt');
    Signal = fscanf(fid, '%g');
    Time=dlmread('TIME.txt');
    fclose(fid);
    window=550*1.38;     %micron
    
%    b = Signal(500:2000); % part of raw data
%    plot(b);

% Digital Filtering:
     Wavelength_Laser=0.6328;     %micron
    Lower_Band=20;
    Upper_Band=135680/2;                %fix(length(Spectrum)/2);
     
     Spectrum=ifft(Signal);
     Spectrum(1:Lower_Band)=0;
     Spectrum(Upper_Band:end)=0;


     Signal_New=fft(Spectrum,[],1);
     Phase_original=angle(Signal_New);
     Phase=unwrap(Phase_original);
    
     
% plot(Phase);     
     Position=-1*Wavelength_Laser*Phase/(2*pi)/2;
     
     dPosition=diff(Position);  %micron
     dTime=diff(Time);          %ms
     
     Velocity_inst=dPosition./dTime;
     Velocity_inst(length(Velocity_inst)+1)=Velocity_inst(length(Velocity_inst));
     
     Velocity_mean=max(Position)/max(Time);
     window_size_time=round(window/Velocity_mean/mean(dTime));
     Velocity_ave=smooth(Velocity_inst,window_size_time);
     plot(Time,Velocity_ave); 
     xlim([Time(window_size_time) Time(length(Time)-window_size_time)]);
     xlabel('Time (ms)');
     ylabel('Stage Velocity (mm/sec)');
     
     Velocity_std=std(Velocity_ave(window_size_time:(length(Velocity_ave)-window_size_time)));
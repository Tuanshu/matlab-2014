clear all

N=1000000;        %入射之光子數(單位時間內, 在特定波長)
IE=0.5;           %interference efficiency

R_sample=0.01;   %assumed
R_stray=R_sample*10;   %assumed

R1_1D=0:0.002:1;        %the reference glass
R2_1D=0:0.002:1;        %the BS


%
R1=repmat(R1_1D,[length(R2_1D) 1]);
R2=repmat(R2_1D',[1 length(R1_1D)]);
T1=1-R1;
T2=1-R2;
% Transmittance

T_Sample=T1.*T2.*R_sample.*T2.*T1;  %the overall T of Mirau for Sample path
T_Reference=T1.*R2.*R1.*R2.*T1;     %the overall T of Mirau for Reference path
T_StrayLight=T1.*T2.*R_stray.*T2.*T1;  %the overall T of Mirau for Sample path
T_Total=T_Sample+T_Reference+T_StrayLight;



N_Sample=T_Sample*N;
N_Reference=T_Reference*N;
N_StrayLight=T_StrayLight*N;


N_Total=N_Sample+N_Reference+N_StrayLight;

N_Interference=2*IE.*(N_Sample.*N_Reference).^0.5;

N_Noise=(N_Total).^0.5;

SNR=(N_Interference./N_Noise).^2;

Equivalent_Interference_Efficiency=N_Interference./N_Total;  %when N_StrayLight=0 and N_Sample=N_Reference, IE=1, this value is equal to 1

imagesc(T_Total,'xdata',R1_1D,'ydata',R2_1D);
xlabel('1st R');
ylabel('BS R');
colormap(gray);
colorbar
axis equal
xlim([0 max(R1_1D)]);
ylim([0 max(R2_1D)]);
caxis([0 0.5]);
%%
imagesc(T_Sample.*T_Reference,'xdata',R1_1D,'ydata',R2_1D);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);
%
imagesc(((T_Sample.*T_Reference).^0.5)./(T_Sample+T_Reference),'xdata',R1_1D,'ydata',R2_1D);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);

%similar to SNR
%imagesc(((T_Sample.*T_Reference))./(T_Sample+T_Reference),'xdata',R1_1D,'ydata',R2_1D);
%xlabel('1st R');
%ylabel('BS R');
%axis equal
%colormap(gray);


imagesc(Equivalent_Interference_Efficiency,'xdata',R1_1D,'ydata',R2_1D);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);
%
imagesc(SNR,'xdata',R1_1D,'ydata',R2_1D);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);


%% Try to consider the CCD saturation by NORMALIZING the system
%簡單說就是:
%1. 先找系統穿透率(T_Sample+T_Reference+T_StrayLight
%2. 決定要在Detector收到多少電子 (N_Detector_Set)
%3. 計算因此對應的入射電子數N_Set (注意! 這個N_Set就會是2D array了)
%4. 用這個N_Set來進行SNR_Set計算

N_Detector_Set=100000;

N_Set=N_Detector_Set./T_Total;
imagesc(N_Set);
caxis([0 100000000]);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);


N_Sample_Set=T_Sample.*N_Set;
N_Reference_Set=T_Reference.*N_Set;
N_StrayLight_Set=T_StrayLight.*N_Set;

N_Total_Set=N_Sample_Set+N_Reference_Set+N_StrayLight_Set;
imagesc(N_Total_Set);
caxis([0 100000000]);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);


N_Interference_Set=2*IE.*(N_Sample_Set.*N_Reference_Set).^0.5;

N_Noise_Set=(N_Total_Set).^0.5;

SNR_Set=(N_Interference_Set./N_Noise_Set).^2;

imagesc(SNR_Set,'xdata',R1_1D,'ydata',R2_1D);
xlabel('1st R');
ylabel('BS R');
axis equal
colormap(gray);


%imagesc(N_Reference_Set./N_Sample_Set,'xdata',R1_1D,'ydata',R2_1D);
%xlabel('1st R');
%ylabel('BS R');
%axis equal
%colormap(gray);
%caxis([0 1000])
%imagesc(N_Set,'xdata',R1_1D,'ydata',R2_1D);
%xlabel('1st R');
%ylabel('BS R');
%axis equal
%colormap(gray);
%caxis([0 N_Detector_Set*100])
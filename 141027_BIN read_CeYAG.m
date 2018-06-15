%% binary to image for FFOCT-1
clear all

ROI=[9 488 9 648];
X_range=[296 405];
Y_range=[201 310];
Z_range=[421 670]; %Ce:YAG glass;%[421 670] Ce:YAG cell;

cd('D:\Users\TuanShu\141022_CeYSO on BCC\bin_CeYAG_125fps_124fpp_0.2micronsecond_log_after realignment_position 3_forget to put LPF\');

N_frame=1971;

% fin=fopen([sprintf('%08d',(1000)),'.bin']);
%  A=(fread(fin,[648 488],'uint16','b'))';
%MAX=0;
%for slice_num=1:N_frame
% fin=fopen([sprintf('%08d',(slice_num)),'.bin']);
% A=(fread(fin,[648 488],'*uint32','b'))';
% B=A(ROI(1):ROI(2),ROI(3):ROI(4));
% if max(max(B))>MAX
%     MAX=max(max(B));
% end
% fclose('all');
% disp(slice_num);
%end

Image_Stack(1:648,1:488,1:N_frame)=0;

for slice_num=1:N_frame
 fin=fopen([sprintf('%08d',(slice_num)),'.bin']);
 Image_Stack(:,:,slice_num)=(fread(fin,[648 488],'float32',0,'b'));
 disp(slice_num);
 fclose('all');
end

Image_Stack_sub=Image_Stack(X_range(1):X_range(2),Y_range(1):Y_range(2),Z_range(1):Z_range(2));
MAX=max(max(max(Image_Stack_sub)));
Image_Stack_sub=Image_Stack_sub/MAX;
%%
N=35;

imagesc(Image_Stack_sub(:,:,N));
colormap(gray);


%%
cd('D:\Users\TuanShu\141020_CeYSO BCC\Bilateral Filtering');

W=10;
SIGMA=[2 0.9];
Image_Stack_BI=Image_Stack_sub;
for p=65:65%1:size(Image_Stack_sub,3)
    temp=Image_Stack_sub(:,:,p);
    Image_Stack_BI(:,:,p)=bfilter2(temp,W,SIGMA);
    disp(p);
end
%
N=65;

imagesc(Image_Stack_BI(:,:,N));
colormap(gray);
axis([equal]);
caxis([0 0.2]);

%%

Scanning_Speed=0.2; %micron/sec
Frame_Rate=125;
FPC=124;

C_per_second=Frame_Rate/FPC;

Micron_per_C=Scanning_Speed/C_per_second;

Position=(1:size(Image_Stack_sub,3))*Micron_per_C;

N=49;
Temp(:,:)=Image_Stack_sub(N,:,:);



for p=1:size(Temp,1)
    [maxvalue maxindex]=max(Temp(p,:));

    Temp_aligned(p,:)=circshift(Temp(p,:),[0 -1*maxindex+50]);
end

imagesc(Temp);
colormap(gray);

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        
Mean_Trace=mean(Temp_aligned(:,:),1)/max(mean(Temp_aligned(:,:),1));

plot(Position-Position(50),Mean_Trace,'-',Position-Position(50),Mean_Trace,'o');
xlabel('Position (\mum)');
ylabel('Normalized Signal Intensity (norm.)');
xlim([-10 10]);


FWHM=Position(find(Mean_Trace>0.5,1,'first'))-Position(find(Mean_Trace>0.5,1,'last'));

%%

interplolation_ratio=10;

Position_interp=interp1(1:length(Position),Position,(1:length(Position)*interplolation_ratio)/interplolation_ratio);
Mean_Trace_interp=interp1(1:length(Mean_Trace),Mean_Trace,(1:length(Mean_Trace)*interplolation_ratio)/interplolation_ratio);

FWHM=abs(Position_interp(find(Mean_Trace_interp>0.5,1,'first'))-Position_interp(find(Mean_Trace_interp>0.5,1,'last')));

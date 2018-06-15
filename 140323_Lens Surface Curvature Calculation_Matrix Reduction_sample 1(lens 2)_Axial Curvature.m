clear all

d_Position=0.56/2/28.75;
Lateral_Spacing=2500/(565-40)/2;

%% Height generation
cd('D:\Users\TuanShu\140304\lens_2\');
Averaging_Factor=1;
Binning_Factor=1;

SPF=50;
LPF=110;   %pixel

image_index=1:2200;
ROI=[1 640;1 512];      %up, down, left, right


Position=(image_index-min(image_index))*d_Position/Averaging_Factor;
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Image_Temp_2(1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
for p=1:length(image_index)
    Image_Temp=dlmread(sprintf(' 140303_les_2_%i',image_index(p)));   
    [m,n]=size(Image_Temp); %M is the original matrix

    Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]),1);
    Image_Temp=reshape(Image_Temp,m/Binning_Factor,[]).';

    Image_Temp=sum(reshape(Image_Temp,Binning_Factor,[]) ,1);
    Image_Temp=reshape(Image_Temp,n/Binning_Factor,[]).';
    if rem(p,Averaging_Factor)==0
        Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2))+Image_Temp_2;
        Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp_2;
        Image_Temp_2=0;
    else
        Image_Temp_2=Image_Temp(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2))+Image_Temp_2;
    end
    %Image_Stack(ceil(p/Averaging_Factor),:,:)=Image_Temp;
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    disp(p);
end

plot(Image_Stack(:,620/2,512/2));


%%
FFT_data=fft(Image_Stack,[],1);
plot(real(FFT_data(:,620/2,512/2)));
%FFT_data(round(length(image_index)/2+1):end,:,:)=0;
FFT_data(LPF:end,:,:)=0;
FFT_data(1:SPF,:,:)=0;
plot(real(FFT_data(:,620/2,512/2)));
data_envelope_temp=abs(ifft(FFT_data));
data_carrier_temp=real(ifft(FFT_data));

plot(Position,data_envelope_temp(:,620/2,512/2),Position,data_carrier_temp(:,620/2,512/2));


%% mask
threshold=25;
[max_value max_index]=max(data_envelope_temp);%data_envelope);
max_intensity(:,:)=max_value;
max_intensity(isnan(max_intensity))=0;
Mask=max_intensity;
Mask(:,:)=1;
Mask(max_intensity<threshold)=0;

imagesc(Mask);

%%
data_envelope=data_envelope_temp;
%% to smooth the 3D data
%smooth_window_3D=1;

%data_envelope=smooth3(data_envelope_temp,'gaussian',smooth_window_3D);

%% 3D resampling
%ratio=1;

%data_envelope_temp=data_envelope;

%xx=1:size(data_envelope_temp,1);
%yy=1:size(data_envelope_temp,2);
%zz=1:size(data_envelope_temp,3);

%xx2=1:ratio:size(data_envelope_temp,1);
%yy2=1:ratio:size(data_envelope_temp,2);
%zz2=1:ratio:size(data_envelope_temp,3);

%data_envelope2=interp3(xx,yy,zz,data_envelope_temp,xx2,yy2,zz2);


%% Take Centriod
%Surface_profile_temp(:,:)=max_index;
clear Surface_profile_temp
index_volume(1:length(image_index),size(data_envelope,2),size(data_envelope,3))=0;
for p=1:length(image_index)
    index_volume(p,:,:)=image_index(p);
end

weighted_data_envelope=data_envelope.*index_volume;
normalized_index=mean(weighted_data_envelope)./mean(data_envelope);
Surface_profile_temp(:,:)=normalized_index;
Surface_profile_temp(isnan(Surface_profile_temp))=0;

imagesc(Surface_profile_temp);
caxis([100 550]);
%% Ave of min and max within a20 window
WindowSize=20;
filter=fspecial('gaussian',WindowSize,round(WindowSize/2));
filter=filter/sum(sum(filter));
Surface_profile_Ave=conv2(Surface_profile_temp,filter,'same');

%% 
Reduction_factor=1;
clear X_temp Y_temp
X_temp(1:round(size(Surface_profile_Ave,1)/Reduction_factor),1:round(size(Surface_profile_Ave,2)/Reduction_factor))=0;
Y_temp(1:round(size(Surface_profile_Ave,1)/Reduction_factor),1:round(size(Surface_profile_Ave,2)/Reduction_factor))=0;

for p=1:size(X_temp,2)
    X_temp(:,p)=p.*Reduction_factor*Lateral_Spacing;
end
for q=1:size(Y_temp,1)
    Y_temp(q,:)=q.*Reduction_factor*Lateral_Spacing;
end

X(1:size(Surface_profile_Ave,1),1:size(Surface_profile_Ave,2))=0;
Y(1:size(Surface_profile_Ave,1),1:size(Surface_profile_Ave,2))=0;

for p=1:size(X,2)
    X(:,p)=p.*Lateral_Spacing;
end
for q=1:size(Y,1)
    Y(q,:)=q.*Lateral_Spacing;
end

Surface_profile_temp=interp2(X,Y,Surface_profile_Ave,X_temp,Y_temp,'nearest');
Surface_profile_temp(isnan(Surface_profile_temp))=0;

%Surface_profile=interp2(X_temp,Y_temp,Surface_profile_temp,X,Y,'spline');

Surface_profile=Surface_profile_temp;

Mask_new=interp2(X,Y,Mask,X_temp,Y_temp,'cubic');
X=X_temp;
Y=Y_temp;
%Mask_new=Mask;
imagesc(Surface_profile);

%%
%plot(1:size(Surface_profile,2),Surface_profile_temp(100,:),1:size(Surface_profile,2),Surface_profile(100,:));
%plot(Surface_profile_reduced(round(100/Reducing_factor),:));

%Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask)*(-1);
Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask_new);
Z(isnan(Z))=0;


imagesc(Z.*Mask,'xdata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,2)*Lateral_Spacing*Reduction_factor),'ydata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,1)*Lateral_Spacing*Reduction_factor));
axis equal
caxis([-12 0]);

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing*Reduction_factor]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing*Reduction_factor]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

%% Axial Curvature Calculation
X_offset=-31;
Y_offset=29;
Z((Mask==0))=min(Z(:));
[maxvalue maxindex_1D]=max(Z(:));
[index_X_Max index_Y_Max]=ind2sub([size(Z,1) size(Z,2)],maxindex_1D);   %use max index as axis

R=((X-X(1,index_X_Max+X_offset)).^2+(Y-Y(index_Y_Max+Y_offset,1)).^2).^0.5;
imagesc(R);

%Sin_Theta=R./(R.^2+Z.^2);
%Theta=asin(Sin_Theta);

Theta=atan(R./abs(Z));
Sin_pi_minus_2Theta=sin(pi-2*Theta);
Axial_Radius_of_Curvature=(R./Sin_pi_minus_2Theta).*Mask/1000;

imagesc(Axial_Radius_of_Curvature,'xdata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,2)*Lateral_Spacing*Reduction_factor),'ydata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,1)*Lateral_Spacing*Reduction_factor));
%plot(Axial_Radius_of_Curvature(300,:))
caxis([10 20]);
axis equal

xlim([Lateral_Spacing*Reduction_factor size(Z,2)*Lateral_Spacing*Reduction_factor]);
ylim([Lateral_Spacing*Reduction_factor size(Z,1)*Lateral_Spacing*Reduction_factor]);
colorbar
xlabel('(micron)');
ylabel('(micron)');
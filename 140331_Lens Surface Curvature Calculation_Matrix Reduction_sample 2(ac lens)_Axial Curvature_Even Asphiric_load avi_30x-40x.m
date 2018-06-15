clear all

d_Position=0.56/2/49;
Lateral_Spacing=122/62;

%% Height generation
Averaging_Factor=1;
Binning_Factor=1;

SPF=35;
LPF=90;   %pixel
%% Height generation
cd('E:\Users\TuanShu\140331_Crystalvue');
MOVIE=aviread('140331_30x-40x_2.avi');
%%
imagesc(MOVIE(191).cdata);
imagetest=MOVIE(191).cdata;
image_index=1:3000;

ROI=[1 344;129 896];      %up, down, left, right


Position=(image_index-min(image_index))*d_Position/Averaging_Factor;
Image_Stack(1:ceil(length(image_index)/Averaging_Factor),1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
Image_Temp_2(1:(ROI(1,2)-ROI(1,1)+1)/Binning_Factor,1:(ROI(2,2)-ROI(2,1)+1)/Binning_Factor)=0;
for p=1:length(image_index)
    Image_Temp=MOVIE(image_index(p)).cdata;
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
plot(Image_Stack(:,344/2,512));

plot(Position,Image_Stack(:,344/2,512));


%%
FFT_data=fft(Image_Stack,[],1);
plot(real(FFT_data(:,344/2,512)));
%FFT_data(round(length(image_index)/2+1):end,:,:)=0;
FFT_data(LPF:end,:,:)=0;
FFT_data(1:SPF,:,:)=0;
plot(real(FFT_data(:,344/2,512)));
data_envelope_temp=abs(ifft(FFT_data));
plot(Position,data_envelope_temp(:,344/2,512));


%% mask
threshold=5;
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
Surface_profile_temp=Surface_profile_temp-max(max(Surface_profile_temp));
imagesc(Surface_profile_temp);
caxis([100 550]);
%% Ave of min and max within a20 window
WindowSize=5;
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
imagesc(Surface_profile.*Mask_new);

%%
%plot(1:size(Surface_profile,2),Surface_profile_temp(100,:),1:size(Surface_profile,2),Surface_profile(100,:));
%plot(Surface_profile_reduced(round(100/Reducing_factor),:));

%Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask)*(-1);
Z=(Surface_profile*(-1).*d_Position.*Mask_new);
Z=Z-max(max(Z));
Z=(Z.*Mask_new);

Z(isnan(Z))=0;

imagesc(Z,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
axis equal
%caxis([0 12]);

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');


%%

%% Axial Curvature Calculation

X_Tilt=-0.0025;
Y_Tilt=0.002;

X_mod(1:size(Surface_profile_Ave,1),1:size(Surface_profile_Ave,2))=0;
Y_mod(1:size(Surface_profile_Ave,1),1:size(Surface_profile_Ave,2))=0;
for p=1:size(X,2)
    X_mod(:,p)=p.*X_Tilt;
end
for q=1:size(Y,1)
    Y_mod(q,:)=q.*Y_Tilt;
end

Z_mod=Z-X_mod-Y_mod;
Z_mod=Z_mod.*Mask_new;
Z_mod=Z_mod-max(max(Z_mod));
imagesc(Z_mod,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
axis equal
%caxis([0 12]);

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

%%


X_offset=430;
Y_offset=198;

R=((X-X(1,X_offset)).^2+(Y-Y(Y_offset,1)).^2).^0.5;
imagesc(R);

%Sin_Theta=R./(R.^2+Z.^2);
%Theta=asin(Sin_Theta);

Theta=atan(R./abs(Z_sum));
Sin_pi_minus_2Theta=sin(pi-2*Theta);
Axial_Radius_of_Curvature=(R./Sin_pi_minus_2Theta).*Mask/1000;

imagesc(R,'xdata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,2)*Lateral_Spacing*Reduction_factor),'ydata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,1)*Lateral_Spacing*Reduction_factor));

Diopter=(1.3375-1)*(1000)./Axial_Radius_of_Curvature;

imagesc(Diopter,'xdata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,2)*Lateral_Spacing*Reduction_factor),'ydata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,1)*Lateral_Spacing*Reduction_factor));
axis equal
xlim([Lateral_Spacing*Reduction_factor size(Z,2)*Lateral_Spacing*Reduction_factor]);
ylim([Lateral_Spacing*Reduction_factor size(Z,1)*Lateral_Spacing*Reduction_factor]);
caxis([0 10]);

colorbar
xlabel('(micron)');
ylabel('(micron)');


imagesc(Axial_Radius_of_Curvature,'xdata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,2)*Lateral_Spacing*Reduction_factor),'ydata',Lateral_Spacing*Reduction_factor:Lateral_Spacing*Reduction_factor:(size(Z,1)*Lateral_Spacing*Reduction_factor));
%plot(Axial_Radius_of_Curvature(300,:))
caxis([0 30]);
axis equal

xlim([Lateral_Spacing*Reduction_factor size(Z,2)*Lateral_Spacing*Reduction_factor]);
ylim([Lateral_Spacing*Reduction_factor size(Z,1)*Lateral_Spacing*Reduction_factor]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

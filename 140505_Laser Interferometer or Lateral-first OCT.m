clear all

Pixel_size=5.5;
M=1;
Lateral_Spacing=Pixel_size*M;

%% Height generation
N=4;
cd('D:\Users\TuanShu\140506_Artificial Eye');
Image=imread(sprintf('%d.png',N))-imread('5.png')-imread('6.png');

imagesc(Image);
colormap(gray);
%%
Y=repmat([1:size(Image,1)]',1,size(Image,2));
X=repmat(1:size(Image,2),size(Image,1),1);



%% Transfer to Polar Cord.
Wavelength=0.635;   %micron
X_Center=600;
Y_Center=920;
R_consider=500;
NNN=2000;
LF=4;

R=repmat([R_consider/NNN:R_consider/NNN:R_consider],NNN,1);
THETA=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);

X_wish=R.*cos(THETA)+X_Center;
Y_wish=R.*sin(THETA)+Y_Center;

Image_Polar_Original=interp2(X,Y,Image,X_wish,Y_wish);
Image_Polar_Original(isnan(Image_Polar_Original))=0;
imagesc(Image_Polar_Original);

%%
WindowSize=200;


filter=fspecial('gaussian',WindowSize,round(WindowSize/2));
for p=1:size(filter,1)
    if p~=round(size(filter,1)/2)
        filter(:,p)=repmat(0,size(filter,2),1);
    end
end
imagesc(filter);

filter=filter/sum(sum(filter));
Image_Polar=conv2(Image_Polar_Original,filter,'same');
imagesc(Image_Polar);
axis equal
%%

Image_Polar_RFFT=fft(Image_Polar,[],2);
imagesc(Image_Polar);
imagesc(real(Image_Polar_RFFT));


plot(real(Image_Polar_RFFT(1000,:)));

Image_Polar_RFFT(:,1:LF)=0;
Image_Polar_RFFT(:,(size(Image_Polar_RFFT,2)/2+1):end)=0;
Image_Polar_Filtered=ifft(Image_Polar_RFFT,[],2);
imagesc(angle(Image_Polar_Filtered));
Phase_Polar=unwrap(angle(Image_Polar_Filtered),[],2);
Phase_Polar(isnan(Phase_Polar))=0;

%% To set equal phase level
R_equal_phase=800;
Phase_Ref=Phase_Polar(size(Phase_Polar,1)/2,R_equal_phase);

Phase_Polar_Equal=Phase_Polar;
for p=1:size(Phase_Polar,1)
    Phase_Polar_Equal(p,:)=Phase_Polar(p,:)+Phase_Ref-Phase_Polar(p,R_equal_phase);
    
end

R_XY=((X-X_Center).^2+(Y-Y_Center).^2).^0.5;
THETA_XY=atan2((Y-Y_Center),(X-X_Center))+pi;

Phase=interp2(R,THETA,Phase_Polar_Equal,R_XY,THETA_XY);

subplot(1,3,1);
imagesc(angle(Image_Polar_Filtered));
axis equal

subplot(1,3,2);
imagesc(Phase);
axis equal

subplot(1,3,3);
imagesc(Image);
axis equal

%%

Z=Phase*Wavelength/2/pi;
Z=Z-min(min(Z));
imagesc(Z);
%%


X_Tilt=0;
Y_Tilt=0;

X_mod(1:size(Z,1),1:size(Z,2))=0;
Y_mod(1:size(Z,1),1:size(Z,2))=0;
for p=1:size(X,2)
    X_mod(:,p)=p.*X_Tilt;
end
for q=1:size(Y,1)
    Y_mod(q,:)=q.*Y_Tilt;
end

Z_mod=Z-X_mod-Y_mod;
%Z_mod=Z_mod.*Mask_new;
Z_mod=Z_mod-max(max(Z_mod));
imagesc(Z_mod,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
axis equal
%caxis([0 12]);

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');
%% Axial Curvature Calculation

R=((X-X(1,X_Center)).^2+(Y-Y(Y_Center,1)).^2).^0.5;
imagesc(R);

%Sin_Theta=R./(R.^2+Z.^2);
%Theta=asin(Sin_Theta);

Theta=atan(R./abs(Z_mod));
Sin_pi_minus_2Theta=sin(pi-2*Theta);
Axial_Radius_of_Curvature=(R./Sin_pi_minus_2Theta)./1000;

imagesc(R,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));

Diopter=(1.3375-1)*(1000)./Axial_Radius_of_Curvature;

imagesc(Diopter,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
axis equal
xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
caxis([0 20]);

colorbar
xlabel('(micron)');
ylabel('(micron)');


imagesc(Axial_Radius_of_Curvature,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
%plot(Axial_Radius_of_Curvature(300,:))
caxis([0 30]);
axis equal

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

Diopter(isinf(Diopter))=14.6;

WindowSize=1;
filter=fspecial('gaussian',WindowSize,round(WindowSize));
filter=filter/sum(sum(filter));
Diopter_Ave=conv2(Diopter,filter,'same');

imagesc(Diopter_Ave,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
%plot(Axial_Radius_of_Curvature(300,:))
caxis([0 200]);
axis equal

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

Diopter_Ave_Mask=Diopter_Ave;
Axial_Radius_of_Curvature_Ave_Mask=(1.3375-1)*(1000)./Diopter_Ave_Mask;


imagesc(Axial_Radius_of_Curvature_Ave_Mask,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
%plot(Axial_Radius_of_Curvature(300,:))
caxis([0 30]);
axis equal

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');
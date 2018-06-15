clear all

d_Position=0.55/2/28.75;
Lateral_Spacing=2500/(579-50)/2;

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
plot(Position,data_envelope_temp(:,620/2,512/2));


%% mask
threshold=50;
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
WindowSize=40;
filter=ones(WindowSize)/WindowSize/WindowSize;
Surface_profile_Ave=conv2(Surface_profile_temp,filter,'same');

%% 
Reduction_factor=8;
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

Surface_profile_temp=interp2(X,Y,Surface_profile_Ave,X_temp,Y_temp,'cubic');
Surface_profile_temp(isnan(Surface_profile_temp))=0;

%Surface_profile=interp2(X_temp,Y_temp,Surface_profile_temp,X,Y,'cubic');

Surface_profile=Surface_profile_temp;

Mask_new=interp2(X,Y,Mask,X_temp,Y_temp,'cubic');
X=X_temp;
Y=Y_temp;
imagesc(Surface_profile);

%%
%plot(1:size(Surface_profile,2),Surface_profile_temp(100,:),1:size(Surface_profile,2),Surface_profile(100,:));
%plot(Surface_profile_reduced(round(100/Reducing_factor),:));

%Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask)*(-1);
Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask_new)*(-1);
Z(isnan(Z))=0;
%imagesc(Z','xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing));

%imagesc(Z)
%axis equal
%xlim([0 size(Z,1)*Lateral_Spacing]);
%ylim([0 size(Z,2)*Lateral_Spacing]);
%xlabel('X Position (mm)');
%ylabel('Y Position (mm)');
%X_1D=X(768/2,:);
%Z_1D=Z(768/2,:);

%plot(X_1D,Z_1D);
%axis equal

%%
%n_cornea=1.3375;

%% Fitting
fitting=0;
if fitting==1
FitSurface= fittype( @(c,r,a, b, x, y) c+(r^2-(x-a).^2-(y-b).^2).^0.5, 'independent', {'x', 'y'},'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2
                                                                                                                        %z=C+(R^2-(x-A)^.2-(y-B)^.2).^0.5
                                                                                                 
Starting_X= 576.1504;                                                                                                                                                                                                                              
Starting_Y=641.1074;
Starting_R=1.0103e+04;
Starting_C=-1.0103e+04;
FitPara=fit([X(:),Y(:)],-Z(:),FitSurface,'Weight',Mask_new(:),'StartPoint',[Starting_C,Starting_R,Starting_X,Starting_Y]);

CheckCurve=FitPara.c+((FitPara.r)^2-(X-FitPara.a).^2-(Y-FitPara.b).^2).^0.5;
CheckCurve(imag(CheckCurve)~=0)=0;

Residual=Z-CheckCurve;

imagesc(CheckCurve);

imagesc(Residual);
caxis([-20 10]);
%%                        
                                                                                                   
Starting_X=500;                                                                                                                                                                                                                              
Starting_Y=600;
Starting_R=10050;
Starting_C=-10040;
imagesc(-Z)
caxis([0 9]);
CheckCurveq=Starting_C+((Starting_R)^2-(X-Starting_X).^2-(Y-Starting_Y).^2).^0.5;
imagesc(real(CheckCurveq));
caxis([0 9]);
end
%%


%X_1D=Lateral_Spacing:Lateral_Spacing:size(Z,2)*Lateral_Spacing;
%Z_1D=Z(100,:)-3.9E-3;

%plot(X_1D,Z_1D);


%Z_Patial=Z(50:size(Z,1)-49,50:size(Z,2)-49);
%surface(-Z_Patial,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Z_Patial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Z_Patial,1)-1));
%surface(-Z,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Z,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Z,1)-1));

%imagesc(-Z,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Z_Patial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Z_Patial,1)-1));
%xlabel('(Micron)');
%ylabel('(Micron)');
%zlabel('(Micron)');

%shading interp
%view(3)

%%

%Z=-CheckCurve;
%Z=Residual;

%%

%[K,H,P1,P2] = surfature(X,Y,Z); 
%X=X;Y=Y;Z=Z;

% First Derivatives
[Xu,Xv] = gradient(X);
[Yu,Yv] = gradient(Y);
[Zu,Zv] = gradient(Z);

% Second Derivatives
[Xuu,Xuv] = gradient(Xu);
[Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu);

[Xuv,Xvv] = gradient(Xv);
[Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv);

% Reshape 2D Arrays into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv(:); 

Xu          =   [Xu Yu Zu];
Xv          =   [Xv Yv Zv];
Xuu         =   [Xuu Yuu Zuu];
Xuv         =   [Xuv Yuv Zuv];
Xvv         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(Xu,Xu,2);
F           =   dot(Xu,Xv,2);
G           =   dot(Xv,Xv,2);

m           =   cross(Xu,Xv,2);
p           =   sqrt(dot(m,m,2));
n           =   m./[p p p]; 

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(Xuu,n,2);
M           =   dot(Xuv,n,2);
N           =   dot(Xvv,n,2);

[s,t] = size(Z);

% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
%K = reshape(K,s,t);

% Mean Curvature
H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
%H = reshape(H,s,t);

% Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);

% eigenvalue: Pmax and Pmin, Matrix

%t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 

%
clear X_P1 Y_P1 Z_P1 X_P2 Y_P2 Z_P2 Curvature_P1 Curvature_P2
X_P1(1:length(Xu))=0;
Y_P1(1:length(Xu))=0;
Z_P1(1:length(Xu))=0;

X_P2(1:length(Xu))=0;
Y_P2(1:length(Xu))=0;
Z_P2(1:length(Xu))=0;

Curvature_P1(1:length(Xu))=0;
Curvature_P2(1:length(Xu))=0;

for p=1:length(Xu)
    
    I=[E(p) F(p); F(p) G(p)];
    II=[L(p) M(p); M(p) N(p)];
    
    MATRIX=I\II;
    [P_dir P_v]=eig(MATRIX);

    %P_dir_1=P_dir(:,1); %the true principal direction is P_dir_1(1)*Xu+P_dir_1(2)*Xv
    %P_dir_2=P_dir(:,2);
    
    Vector_1=[P_dir(1,1)*Xu(p,1)+P_dir(2,1)*Xv(p,1) P_dir(1,1)*Xu(p,2)+P_dir(2,1)*Xv(p,2) P_dir(1,1)*Xu(p,3)+P_dir(2,1)*Xv(p,3)];
    Vector_2=[P_dir(1,2)*Xu(p)+P_dir(2,2)*Xv(p,1) P_dir(1,2)*Xu(p,2)+P_dir(2,2)*Xv(p,2) P_dir(1,2)*Xu(p,3)+P_dir(2,2)*Xv(p,3)];
    Vector_1=Vector_1/norm(Vector_1);
    Vector_2=Vector_2/norm(Vector_2);
    
    X_P1(p)=Vector_1(1);
    Y_P1(p)=Vector_1(2);
    Z_P1(p)=Vector_1(3);
    
    X_P2(p)=Vector_2(1);
    Y_P2(p)=Vector_2(2);
    Z_P2(p)=Vector_2(3);
    %X_P1(p)=P_dir(1,1)*Xu(p)+P_dir(2,1)*Xv(p);
    %Y_P1(p)=P_dir(1,1)*Yu(p)+P_dir(2,1)*Yv(p);
    %Z_P1(p)=P_dir(1,1)*Zu(p)+P_dir(2,1)*Zv(p);

    %X_P2(p)=P_dir(1,2)*Xu(p)+P_dir(2,2)*Xv(p);
    %Y_P2(p)=P_dir(1,2)*Yu(p)+P_dir(2,2)*Yv(p);
    %Z_P2(p)=P_dir(1,2)*Zu(p)+P_dir(2,2)*Zv(p);
    
    Curvature_P1(p)=-1*P_v(1,1);
    Curvature_P2(p)=-1*P_v(2,2);
    
    disp(p);
    
end

X_P1 = reshape(X_P1,s,t);
Y_P1 = reshape(Y_P1,s,t);
Z_P1 = reshape(Z_P1,s,t);

X_P2 = reshape(X_P2,s,t);
Y_P2 = reshape(Y_P2,s,t);
Z_P2 = reshape(Z_P2,s,t);

Curvature_P1 = reshape(Curvature_P1,s,t);
Curvature_P2 = reshape(Curvature_P2,s,t);
Curvature_P1(isnan(Curvature_P1))=0;
Curvature_P2(isnan(Curvature_P2))=0;

D=Curvature_P1>Curvature_P2;

Pmax = reshape(Pmax,s,t);
Pmin = reshape(Pmin,s,t);

K = reshape(K,s,t);
H = reshape(H,s,t);
Curvature_P1=Curvature_P1.*Mask_new;
Curvature_P2=Curvature_P2.*Mask_new;

imagesc(Curvature_P1);
axis equal
xlabel('Pixel Number');
ylabel('Pixel Number');
caxis([-0.0002 0]);
%%
WindowSize2=10;
filter2=ones(WindowSize2)/WindowSize2/WindowSize2;
Curvature_P1_smooth=conv2(Curvature_P1,filter2,'same').*Mask_new;
Curvature_P2_smooth=conv2(Curvature_P2,filter2,'same').*Mask_new;

R_P1_smooth=-1./(Curvature_P1_smooth)/1000;    %mm
R_P2_smooth=-1./(Curvature_P2_smooth)/1000;    %mm

imagesc(Curvature_P1_smooth);
axis equal
caxis([-0.0002 0]);


imagesc(R_P1_smooth);
caxis([10 20]);
axis equal

%%

Reduction_factor2=1;
clear X_temp2 Y_temp2 X2 Y2
X_temp2(1:round(size(Curvature_P1_smooth,1)/Reduction_factor2),1:round(size(Curvature_P1_smooth,2)/Reduction_factor2))=0;
Y_temp2(1:round(size(Curvature_P1_smooth,1)/Reduction_factor2),1:round(size(Curvature_P1_smooth,2)/Reduction_factor2))=0;

for p=1:size(X_temp2,2)
    X_temp2(:,p)=p*Reduction_factor2;
end
for q=1:size(Y_temp2,1)
    Y_temp2(q,:)=q*Reduction_factor2;
end

X2(1:size(Curvature_P1_smooth,1),1:size(Curvature_P1_smooth,2))=0;
Y2(1:size(Curvature_P1_smooth,1),1:size(Curvature_P1_smooth,2))=0;

for p=1:size(X2,2)
    X2(:,p)=p;
end
for q=1:size(Y2,1)
    Y2(q,:)=q;
end

Curvature_P1_interp_temp=interp2(X2,Y2,Curvature_P1_smooth,X_temp2,Y_temp2,'spline');
Curvature_P1_interp_temp(isnan(Curvature_P1_interp_temp))=0;
Curvature_P1_interp=interp2(X_temp2,Y_temp2,Curvature_P1_interp_temp,X2,Y2,'spline');


R_P1_interp=-1./(Curvature_P1_interp)/1000;    %mm
R_P2_interp=-1./(Curvature_P1_interp)/1000;    %mm

imagesc(R_P1_interp);
caxis([13 16]);
axis equal

%%
axis equal
xlim([0 1024]);
ylim([0 768]);
xlabel('Pixel Number');
ylabel('Pixel Number');
caxis([-10 10]);

%imagesc(Curvature_P1_smooth);
%axis equal
%xlim([0 1024]);
%ylim([0 768]);
%xlabel('Pixel Number');
%ylabel('Pixel Number');
%caxis([-0.01 0.01]);

Downsample_Factor=5;

X_resam=X(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Y_resam=Y(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Z_resam=Z(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);

X_P1_resam=X_P1(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Y_P1_resam=Y_P1(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Z_P1_resam=Z_P1(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);


X_P2_resam=X_P2(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Y_P2_resam=Y_P2(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Z_P2_resam=Z_P2(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);

test=X_P1_resam.*X_P2_resam+Y_P1_resam.*Y_P2_resam+Z_P1_resam.*Z_P2_resam;

quiver3(X_resam,Y_resam,Z_resam,X_P1_resam,Y_P1_resam,Z_P1_resam);
AA=500;
h=quiver3(X_resam,Y_resam,Z_resam,X_P1_resam,Y_P1_resam,Z_P1_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
%h=quiver3(X_resam,Y_resam,Z_resam,Xr_XY_resam,Yr_XY_resam,Zr_XY_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
%adjust_quiver_arrowhead_size(h,2);
Z_graph=Z;
Z_graph(Z_graph<0.007)=NaN;
surface(X,Y,Z_graph,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
surface(X_resam,Y_resam,Z_resam);
xlabel('(Micron)');
ylabel('(Micron)');
zlabel('(Micron)');
shading interp
view(3)
%view(2)
%axis equal
xlim([0 max(max(X))]);
ylim([0 max(max(Y))]);
zlim([0.007 max(max(Z))]);




%R1=1./real(P1)/1000;
%R2=1./real(P2)/1000;
%RK=1./real(K)/1000/1000;
%RH=1./real(H)/1000;



%imagesc(-R1,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(R1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(R1,1)-1));
%xlabel('(Micron)');
%ylabel('(Micron)');
%colorbar
%axis equal
%axis([0 Lateral_Spacing*(size(R1,2)-1) 0 Lateral_Spacing*(size(R1,1)-1)])
%caxis([5 15])
%caxis([0 20])


%% Next I need the radial direction (pre-experimental knowledge)
% Radial direction: in fact I want to find (polar coordinate) the tangent
% vector along R
NNN=2000;
R_consider=Lateral_Spacing*max(size(Z,1),size(Z,2));
Center=[size(X,1)/2 size(X,2)/2]*Lateral_Spacing;


R=repmat([R_consider/NNN:R_consider/NNN:R_consider],NNN,1);
THETA=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);

[Rr Rtheta]=gradient(R);
[THETAr THETAtheta]=gradient(THETA);

X_wish=R.*cos(THETA)+Center(2);
Y_wish=R.*sin(THETA)+Center(1);

X_Polar=interp2(X,Y,X,X_wish,Y_wish);
Y_Polar=interp2(X,Y,Y,X_wish,Y_wish);
Z_Polar=interp2(X,Y,Z,X_wish,Y_wish);

[Xr Xtheta]=gradient(X_Polar);
[Yr Ytheta]=gradient(Y_Polar);
[Zr Ztheta]=gradient(Z_Polar);

% 現在來求
% back to XY coordinate

Center_2=[size(X,1)/2 size(X,2)/2]*Lateral_Spacing;


R_XY=((X-Center_2(2)).^2+(Y-Center_2(1)).^2).^0.5;
THETA_XY=atan2((Y-Center_2(1)),(X-Center_2(2)))+pi;

Xr_XY=interp2(R,THETA,Xr,R_XY,THETA_XY);
Yr_XY=interp2(R,THETA,Yr,R_XY,THETA_XY);
Zr_XY=interp2(R,THETA,Zr,R_XY,THETA_XY);

NORMr_XY=(Xr_XY.^2+Yr_XY.^2+Zr_XY.^2).^0.5;
Xr_XY=Xr_XY./NORMr_XY;
Yr_XY=Yr_XY./NORMr_XY;
Zr_XY=Zr_XY./NORMr_XY.*(-1);      %不是很確定原因, 但這裡Z的分量必須乘個負號, 不然畫出來的vector方向會錯(也會影響到結果?)

%這就是我要的radial方向的tangent vector


Xr_XY_resam=Xr_XY(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Yr_XY_resam=Yr_XY(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
Zr_XY_resam=Zr_XY(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);



%現在來求tangent vector和P1 vector的夾角phi

PHI=acos(X_P1.*Xr_XY+Y_P1.*Yr_XY+Z_P1.*Zr_XY);

Curvature_Radial=(Curvature_P1.*cos(PHI).^2+Curvature_P2.*sin(PHI).^2);
Curvature_Radial(round(size(Curvature_Radial,1)/2)+1,:)=(Curvature_Radial(round(size(Curvature_Radial,1)/2)+2,:)+Curvature_Radial(round(size(Curvature_Radial,1)/2),:))/2;
Curvature_Radial(:,round(size(Curvature_Radial,1)/2)+1)=(Curvature_Radial(:,round(size(Curvature_Radial,1)/2)+2)+Curvature_Radial(:,round(size(Curvature_Radial,1)/2)))/2;
Curvature_Radial(isnan(Curvature_Radial))=0;

imagesc(real(PHI));
axis equal

Curvature_Radial_TEST=(Curvature_P1.*sin(PHI).^2+Curvature_P2.*cos(PHI).^2);
Curvature_Radial_TEST(round(size(Curvature_Radial_TEST,1)/2)+1,:)=(Curvature_Radial_TEST(round(size(Curvature_Radial_TEST,1)/2)+2,:)+Curvature_Radial_TEST(round(size(Curvature_Radial_TEST,1)/2),:))/2;
Curvature_Radial_TEST(:,round(size(Curvature_Radial_TEST,1)/2)+1)=(Curvature_Radial_TEST(:,round(size(Curvature_Radial_TEST,1)/2)+2)+Curvature_Radial_TEST(:,round(size(Curvature_Radial_TEST,1)/2)))/2;
Curvature_Radial_TEST(isnan(Curvature_Radial_TEST))=0;

%imagesc(Curvature_Radial);
%Curvature_center=Curvature_Radial(round(size(Curvature_Radial,1)/2),round(size(Curvature_Radial,2)/2));
%caxis([Curvature_center/2 Curvature_center*1.5]);
%axis equal

Diopter_Radial=(n_cornea-1).*1E6.*Curvature_Radial;
Diopter_Radial_TEST=(n_cornea-1).*1E6.*Curvature_Radial_TEST;

Diopter_Principal_1=(n_cornea-1).*1E6.*Curvature_P1;
Diopter_Principal_2=(n_cornea-1).*1E6.*Curvature_P2;

%subplot(1,2,1)
imagesc(Diopter_Radial,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,1)-1));
caxis([35 45]);
xlabel('(Micron)');
ylabel('(Micron)');
colorbar 
axis equal
xlim([0 8700]);
ylim([0 8700]);

%subplot(1,2,2)
%imagesc(Diopter_Radial_TEST,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial_TEST,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial_TEST,1)-1));
%caxis([30 45]);
%xlabel('(Micron)');
%ylabel('(Micron)');
%colorbar 
%axis equal
%xlim([0 8700]);
%ylim([0 8700]);

%%
principal_plot=0;
if principal_plot==1
    subplot(2,2,1)
    
Diopter_Principal_1(round(size(Diopter_Principal_1,1)/2)+1,:)=(Diopter_Principal_1(round(size(Diopter_Principal_1,1)/2)+2,:)+Diopter_Principal_1(round(size(Diopter_Principal_1,1)/2),:))/2;
Diopter_Principal_1(:,round(size(Diopter_Principal_1,1)/2)+1)=(Diopter_Principal_1(:,round(size(Diopter_Principal_1,1)/2)+2)+Diopter_Principal_1(:,round(size(Diopter_Principal_1,1)/2)))/2;
    
    
    imagesc(Diopter_Principal_1,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,1)-1));
    caxis([35 45]);
    xlabel('(Micron)');
    ylabel('(Micron)');
    colorbar 
    axis equal
    xlim([0 8700]);
    ylim([0 8700]);
    
Diopter_Principal_2(round(size(Diopter_Principal_2,1)/2)+1,:)=(Diopter_Principal_2(round(size(Diopter_Principal_2,1)/2)+2,:)+Diopter_Principal_2(round(size(Diopter_Principal_2,1)/2),:))/2;
Diopter_Principal_2(:,round(size(Diopter_Principal_2,1)/2)+1)=(Diopter_Principal_2(:,round(size(Diopter_Principal_2,1)/2)+2)+Diopter_Principal_2(:,round(size(Diopter_Principal_2,1)/2)))/2;
    
    
    imagesc(Diopter_Principal_2,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,1)-1));
    caxis([35 45]);
    xlabel('(Micron)');
    ylabel('(Micron)');
    colorbar 
    axis equal
    xlim([0 8700]);
    ylim([0 8700]);

    subplot(2,2,2)
    h=quiver3(X_resam,Y_resam,Z_resam,X_P1_resam,Y_P1_resam,Z_P1_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
    %h=quiver3(X_resam,Y_resam,Z_resam,Xr_XY_resam,Yr_XY_resam,Zr_XY_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
    %adjust_quiver_arrowhead_size(h,2);
    surface(X,Y,Z,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
    xlim([0 8700]);
    ylim([0 8700]);
    xlabel('(Micron)');
    ylabel('(Micron)');
    zlabel('(Micron)');
    shading interp
    %view(3)
    view(2)
    axis equal
    
    subplot(2,2,3)
    imagesc(Diopter_Principal_2,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Diopter_Radial,1)-1));
    caxis([20 60]);
    xlabel('(Micron)');
    ylabel('(Micron)');
    colorbar 
    axis equal
    xlim([0 8700]);
    ylim([0 8700]);

    subplot(2,2,4)
    h=quiver3(X_resam,Y_resam,Z_resam,X_P2_resam,Y_P2_resam,Z_P2_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
    %h=quiver3(X_resam,Y_resam,Z_resam,Xr_XY_resam,Yr_XY_resam,Zr_XY_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
    %adjust_quiver_arrowhead_size(h,2);
    surface(X,Y,Z,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
    xlim([0 8700]);
    ylim([0 8700]);
    xlabel('(Micron)');
    ylabel('(Micron)');
    zlabel('(Micron)');
    shading interp
    %view(3)
    view(2)
    axis equal
end
%imagesc(Diopter_Principal_2);
%Curvature_center=Diopter_Radial(round(size(Curvature_Radial,1)/2),round(size(Curvature_Radial,2)/2));
%caxis([35 45]);
%axis equal


%X_Radial=X-X_Center;
%Y_Radial=Y-Y_Center;
%Norm_Radial=(X_Radial.^2+Y_Radial.^2).^0.5;
%X_Radial=X_Radial./Norm_Radial;
%Y_Radial=Y_Radial./Norm_Radial;
%Z?

%X_Radial_resam=X_Radial(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);
%Y_Radial_resam=Y_Radial(Downsample_Factor:Downsample_Factor:end,Downsample_Factor:Downsample_Factor:end);

%quiver3(X_resam,Y_resam,Z_resam,X_P2_resam,Y_P2_resam,Z_P2_resam);


%% trans to polar
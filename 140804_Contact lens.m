clear all

Pixel_size=4.65;
M=1.09;
Lateral_Spacing=Pixel_size*M;
KK=350;
%ROI=[1 344 KK KK+343];
% Height generation
Sample_ave=10;
Sample_strat_index=628;%598;
Background_ave=1;
Background_strat_index=800;
%%
cd('D:\Users\TuanShu\140804_Contact Lens\2 interfaces\');

N1=471;
N2=33;  %%91
dN=1;
Image_test=dlmread(sprintf('interface_1_for_diff__%d',N2))-dlmread(sprintf('interface_1_for_diff__%d',N2+dN));
imagesc(Image_test);
xlim([500 1600]);
ylim([500 1600]);
%%
clear Image_temp Image_background_temp
for N=1:1
cd('D:\Users\TuanShu\140804_Contact Lens\2 interfaces\');
%Image=imread(sprintf('%d.png',N));%imread(sprintf('%d.png',N))-imread('5.png')-imread('6.png');
for p=1:Sample_ave
    Image_temp(:,:,p)=dlmread(sprintf('interface_1_480_%d',p+Sample_strat_index));
    disp(p);
end
for p=1:Background_ave
    Image_background_temp(:,:,p)=dlmread(sprintf('interface_2_464_%d',p+Background_strat_index));
    disp(p);
end
%Image=dlmread('Lower_test_1_160')-dlmread('Lower_test_background_174');%imread(sprintf('%d.png',N))-imread('5.png')-imread('6.png');
Image=dlmread(sprintf('interface_1_for_diff__%d',N2))-dlmread(sprintf('interface_1_for_diff__%d',N2+dN));%mean(Image_temp,3)-mean(Image_background_temp,3);%imread(sprintf('%d.png',N))-imread('5.png')-imread('6.png');
%Image=double(Image(ROI(1):ROI(2),ROI(3):ROI(4)));
imagesc(Image);
colormap(gray);
axis equal
%%
Y=repmat([1:size(Image,1)]',1,size(Image,2));
X=repmat(1:size(Image,2),size(Image,1),1);


%% Transfer to Polar Cord.
Wavelength=0.830;   %micron
X_Center=1106;
Y_Center=1051;
R_consider=850;
NNN=2000;


R=repmat([R_consider/NNN:R_consider/NNN:R_consider],NNN,1);
THETA=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);

X_wish=R.*cos(THETA)+X_Center;
Y_wish=R.*sin(THETA)+Y_Center;

Image_Polar_Original=interp2(X,Y,Image,X_wish,Y_wish);
Image_Polar_Original(isnan(Image_Polar_Original))=0;
imagesc(Image_Polar_Original);


%
LF=6;
HF=135;

Image_Polar_RFFT=fft(Image_Polar_Original,[],2);
imagesc(real(Image_Polar_RFFT));


plot(real(Image_Polar_RFFT(1000,:)));

Image_Polar_RFFT(:,1:LF)=0;
Image_Polar_RFFT(:,HF:end)=0;
%Image_Polar_RFFT(:,(size(Image_Polar_RFFT,2)/2+1):end)=0;
Image_Polar_Filtered=ifft(Image_Polar_RFFT,[],2);
%Image_Polar_Filtered_Envelope=abs(ifft(Image_Polar_RFFT,[],2));

%imagesc(real(Image_Polar_Filtered_Envelope));

%
WindowSize=10;


filter=fspecial('gaussian',WindowSize,round(WindowSize/2));
for p=1:size(filter,1)
    if p~=round(size(filter,1)/2)
        filter(:,p)=repmat(0,size(filter,2),1);
    end
end
imagesc(filter);

filter=filter/sum(sum(filter));
%Image_Polar=conv2(Image_Polar_Original,filter,'same');
Image_Polar=conv2(real(Image_Polar_Filtered),filter,'same');

Image_Polar=max(max(Image_Polar))-Image_Polar;



imagesc(Image_Polar);
axis equal

% select特定等高線
Rough_R=90;%90;%280+200*(N-1);
first_tolerance=10;
Jump=2;
tolerance=1;
number_of_line=15;
Next_Start_Index=1;
clear Max_line_index
for q=1:number_of_line
    for p=1:NNN/Jump
        if (q==1)&&(p==1)
            if Jump==1
                [max_value max_index]=max(Image_Polar(p,(Rough_R-first_tolerance):(Rough_R+first_tolerance)));
                Max_line_index(p,q)=max_index+Rough_R-first_tolerance-1;
            else
                [max_value max_index]=max(Image_Polar((p*Jump),(Rough_R-first_tolerance):(Rough_R+first_tolerance)));
                Max_line_index(((p-1)*Jump+1):(p*Jump),q)=interp1([((p-1)*Jump+1) (p*Jump)],[0 max_index],[((p-1)*Jump+1):(p*Jump)])+Rough_R-first_tolerance-1;
            end
        elseif (p==1)
            if Jump==1
            	[max_value max_index]=max(Image_Polar(p,(Next_Start_Index-tolerance):(Next_Start_Index+tolerance)));
            	Max_line_index(p,q)=max_index+Next_Start_Index-tolerance-1;
            else
                [max_value max_index]=max(Image_Polar((p*Jump),(Next_Start_Index-tolerance):(Next_Start_Index+tolerance)));
                Max_line_index(((p-1)*Jump+1):(p*Jump),q)=interp1([((p-1)*Jump+1) (p*Jump)],[0 max_index],[((p-1)*Jump+1):(p*Jump)])+Next_Start_Index-tolerance-1;
            end
        else
            if Jump==1
                %[max_value max_index]=max(Image_Polar(p,(Max_line_index(p-1,q)-tolerance):(Max_line_index(p-1,q)+tolerance)));
                [max_value max_index]=max(Image_Polar(p,(Max_line_index(p-1,q)-tolerance):(Max_line_index(p-1,q)+tolerance)));
                Max_line_index(p,q)=max_index+Max_line_index(p-1,q)-tolerance-1;
            else
                %[max_value max_index]=max(Image_Polar((p*Jump),(Max_line_index(((p-1)*Jump),q)-tolerance):(Max_line_index(((p-1)*Jump),q)+tolerance)));
                %Max_line_index(((p-1)*Jump):(p*Jump),q)=interp1([((p-1)*Jump) (p*Jump)],[0 max_index],[((p-1)*Jump):(p*Jump)])+Max_line_index(((p-1)*Jump),q)-tolerance-1;

                Start_index=Max_line_index(((p-1)*Jump),q);

                [max_value max_index]=max(Image_Polar((p*Jump),(Start_index-tolerance):(Start_index+tolerance)));
                Max_line_index(((p-1)*Jump):(p*Jump),q)=interp1([((p-1)*Jump) (p*Jump)],[0 max_index],[((p-1)*Jump):(p*Jump)])+Start_index-tolerance-1;
            end
        end
    end
    [next_min_value next_min_index]=findpeaks(max(max(Image_Polar))-Image_Polar(1,(Max_line_index(1,q)):end),'NPEAKS',1);
    [next_man_value next_max_index]=findpeaks(Image_Polar(1,(Max_line_index(1,q)+next_min_index):end),'NPEAKS',1);
    Next_Start_Index=next_min_index+next_max_index+Max_line_index(1,q)-1;
    
end
    %findpeaks(Image_Polar(1,Max_line_index(1,q):end),'NPEAKS',1)+Max_line_index(1,q)-1;

imagesc(Image_Polar);
%axis equal
colormap(gray);

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
hold on
for q=1:number_of_line
    plot(Max_line_index(:,q),1:NNN);
end
hold off

%% back to the original RTheta base
basic_relative_height_array=-1*Wavelength/2*(1:number_of_line);
Z_RTheta(1:NNN,1:NNN)=0;
for p=1:NNN
    Z_RTheta(p,:)=interp1(Max_line_index(p,:),basic_relative_height_array,1:NNN);
end
%Z_RTheta(isnan(Z_RTheta))=0;


R_XY=((X-X_Center).^2+(Y-Y_Center).^2).^0.5;
THETA_XY=atan2((Y-Y_Center),(X-X_Center))+pi;

Z_XY(:,:,N)=interp2(R,THETA,Z_RTheta,R_XY,THETA_XY);
imagesc(Z_XY(:,:,N));
axis equal
colormap(gray);

%clear Max_line_index
end



TFmap_1=not(isnan(Z_XY(:,:,1)));
TFmap_2=not(isnan(Z_XY(:,:,2)));
TFmap_3=not(isnan(Z_XY(:,:,3)));

%%
NN=3;
Z_temp=Z_XY(:,:,NN);
Z_temp(isnan(Z_temp))=0;
imagesc(Z_temp,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z_temp,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z_temp,1)*Lateral_Spacing));
axis equal
xlim([0 size(Z_temp,2)*Lateral_Spacing]);
ylim([0 size(Z_temp,2)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');
colormap(gray);

%%
Sample_index=[60 170];

for N=2:3
    
    Z_XY(:,:,N)=Z_XY(Sample_index(1),Sample_index(2),1)-Z_XY(Sample_index(1),Sample_index(2),N)+Z_XY(:,:,N);

end

plot((1:344)*Lateral_Spacing,Z_XY(1:344,170,1),(1:344)*Lateral_Spacing,Z_XY(1:344,170,2),(1:344)*Lateral_Spacing,Z_XY(1:344,170,3),'LineWidth',2);

xlabel('(micron)');
ylabel('(micron)');
%
Z_XY(isnan(Z_XY))=0;
Z_allsumup=sum(Z_XY,3);
imagesc(Z_allsumup);

TFmap_total=TFmap_1+TFmap_2+TFmap_3;
imagesc(TFmap_total);
%%
Z_allsumup(TFmap_total==2)=Z_allsumup(TFmap_total==2)/2;
Z_allsumup(TFmap_total==3)=Z_allsumup(TFmap_total==3)/3;

imagesc(Z_allsumup,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z_temp,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z_temp,1)*Lateral_Spacing));
axis equal
xlim([0 size(Z_temp,2)*Lateral_Spacing]);
ylim([0 size(Z_temp,2)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');
colormap(gray);
%% Extrapolation
Mask_Z_allsumup=Z_allsumup;
Mask_Z_allsumup(Mask_Z_allsumup~=0)=1;

%%
Z=-1*Z_allsumup;
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
%Z_mod=Z_mod-max(max(Z_mod));
imagesc(Z_mod,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
axis equal
%caxis([0 12]);

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

%% Z_ref

R_Reference=10.3;%9;%7.79;          %mm
conic_Ref=0;

R=Lateral_Spacing*(((X-X(1,X_Center)).^2+(Y-Y(Y_Center,1)).^2).^0.5);
imagesc(R);


Z_ref=((R.^2)/(R_Reference*1000))./(1+(1-(1+conic_Ref).*(R/(R_Reference*1000)).^2).^0.5);

Z_ref(Z_mod==0)=0;

%% Z_total
Z_offset=0;
Z_total=Z_ref+Z_mod;
imagesc(Z_total,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));

%% Axial Curvature
%Sin_Theta=R./(R.^2+Z.^2);
%Theta=asin(Sin_Theta);

Theta=atan(R./abs(Z_total+Z_offset));
Sin_pi_minus_2Theta=sin(pi-2*Theta);
Axial_Radius_of_Curvature=(R./Sin_pi_minus_2Theta)./1000;

imagesc(R,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));

Diopter=(1.3375-1)*(1000)./Axial_Radius_of_Curvature;

imagesc(Diopter,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
axis equal
xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
caxis([40 50]);

colorbar
xlabel('(micron)');
ylabel('(micron)');


imagesc(Axial_Radius_of_Curvature,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
%plot(Axial_Radius_of_Curvature(300,:))
caxis([7 8.5]);
axis equal

xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
colorbar
xlabel('(micron)');
ylabel('(micron)');

%% Principle Curvature

Y_resam=Lateral_Spacing*2*repmat([1:size(Image,1)/2]',1,size(Image,2)/2);
X_resam=Lateral_Spacing*2*repmat(1:size(Image,2)/2,size(Image,1)/2,1);

Z_resam=interp2(Lateral_Spacing*X,Lateral_Spacing*Y,Z,X_resam,Y_resam);
Z_resam(isnan(Z_resam))=0;



Z_final=Z_resam;
X_final=X_resam;
Y_final=Y_resam;
%%

[Xu,Xv] = gradient(X_final);
[Yu,Yv] = gradient(Y_final);
[Zu,Zv] = gradient(Z_final);

% Second Derivatives
[Xuu,Xuv] = gradient(Xu);
[Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu);

[Xuv,Xvv] = gradient(Xv);
[Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv);

%% Surface reconstruction
Smooth_window=10;
Start_pixel_X=104;
Start_pixel_Y=111;

filter=fspecial('gaussian',Smooth_window,round(Smooth_window/2));

Zuu_filtered=conv2(Zuu,filter,'same');
Zuv_filtered=conv2(Zuv,filter,'same');
Zvv_filtered=conv2(Zvv,filter,'same');

Zv_rec(1:size(Zv,1),1:size(Zv,2))=0;
Zu_rec(1:size(Zu,1),1:size(Zu,2))=0;

for p=1:size(Zv,1)
    for q=1:size(Zv,2)
        %Zu_rec(p,q)=Zu_rec(Start_pixel_X,Start_pixel_Y)+sum(Zuv_filtered(Start_pixel_X:p,q),1)+sum(Zuu_filtered(p,Start_pixel_Y:q),2);
        %Zv_rec(p,q)=Zv_rec(Start_pixel_X,Start_pixel_Y)+sum(Zvv_filtered(Start_pixel_X:p,q),1)+sum(Zuv_filtered(p,Start_pixel_Y:q),2);
        if p>Start_pixel_X
            Zu_x_increasement=sum(Zuv_filtered(Start_pixel_X:p,q),1);
            Zv_x_increasement=sum(Zvv_filtered(Start_pixel_X:p,q),1);
        else
            Zu_x_increasement=-1*sum(Zuv_filtered(p:Start_pixel_X,q),1);
            Zv_x_increasement=-1*sum(Zvv_filtered(p:Start_pixel_X,q),1);
        end
        if q>Start_pixel_Y
            Zu_y_increasement=sum(Zuu_filtered(p,Start_pixel_Y:q),2);
            Zv_y_increasement=sum(Zuv_filtered(p,Start_pixel_Y:q),2);
        else
            Zu_y_increasement=-1*sum(Zuu_filtered(p,q:Start_pixel_Y),2);
            Zv_y_increasement=-1*sum(Zuv_filtered(p,q:Start_pixel_Y),2);
        end 
        Zu_rec(p,q)=Zu(Start_pixel_X,Start_pixel_Y)+Zu_x_increasement+Zu_y_increasement;
        Zv_rec(p,q)=Zv(Start_pixel_X,Start_pixel_Y)+Zv_x_increasement+Zv_y_increasement;
    end
end

subplot(2,2,1)
imagesc(Zuu_filtered>0);
caxis([-1 1]);

subplot(2,2,2)
imagesc(Zvv_filtered>0);
caxis([-1 1]);

subplot(2,2,3)
imagesc(Zu);
caxis([-1 1]);

subplot(2,2,4)
imagesc(Zu_rec);
caxis([-1 1]);

%

Z_rec(1:size(Z_final,1),1:size(Z_final,2))=0;

for p=1:size(Zv,1)
    for q=1:size(Zv,2)
        %Zu_rec(p,q)=Zu_rec(Start_pixel_X,Start_pixel_Y)+sum(Zuv_filtered(Start_pixel_X:p,q),1)+sum(Zuu_filtered(p,Start_pixel_Y:q),2);
        %Zv_rec(p,q)=Zv_rec(Start_pixel_X,Start_pixel_Y)+sum(Zvv_filtered(Start_pixel_X:p,q),1)+sum(Zuv_filtered(p,Start_pixel_Y:q),2);
        if p>Start_pixel_X
            Z_x_increasement=sum(Zv_rec(Start_pixel_X:p,q),1);
        else
            Z_x_increasement=-1*sum(Zv_rec(p:Start_pixel_X,q),1);
        end
        if q>Start_pixel_Y
            Z_y_increasement=sum(Zu_rec(p,Start_pixel_Y:q),2);
        else
            Z_y_increasement=-1*sum(Zu_rec(p,q:Start_pixel_Y),2);
        end 
        Z_rec(p,q)=Z_final(Start_pixel_X,Start_pixel_Y)+Z_x_increasement+Z_y_increasement;
        Z_rec(p,q)=Z_final(Start_pixel_X,Start_pixel_Y)+Z_x_increasement+Z_y_increasement;
    end
end

subplot(2,2,1)
imagesc(Z);
caxis([0 10]);

subplot(2,2,2)
imagesc(Z_rec);
caxis([0 10]);

subplot(2,2,3)
imagesc(Zu);
caxis([-1 1]);

subplot(2,2,4)
imagesc(Zu_rec);
caxis([-1 1]);


%%
% Reshape 2D Arrays into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu_rec(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv_rec(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu_filtered(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv_filtered(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv_filtered(:); 

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

[s,t] = size(Z_final);

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
%Curvature_P1=Curvature_P1.*Mask_new;
%Curvature_P2=Curvature_P2.*Mask_new;

Curvature_max=D.*Curvature_P1+(1-D).*Curvature_P2;
Curvature_min=(D-1).*Curvature_P1+D.*Curvature_P2;


imagesc(Curvature_P1);
axis equal
xlabel('Pixel Number');
ylabel('Pixel Number');
caxis([-0.00005 0.00005]);
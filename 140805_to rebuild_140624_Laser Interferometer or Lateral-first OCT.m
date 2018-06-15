clear all

Pixel_size=5.5;
M=0.6;
Lateral_Spacing=Pixel_size*M;
KK=350;
N_Start=1;
N_End=13;
% Height generation

cd('D:\Users\TuanShu\140624\');
%% background

Image_Stack(1:2040,1:2048,1:(N_End-N_Start+1))=0;

for xx=1:(N_End-N_Start+1)

Image_Stack(:,:,xx)=double(dlmread(sprintf('140624_%d',N_Start+xx-1)));
disp(xx);
end


%%

NN=10;

QQ(:,:)=Image_Stack(:,:,NN);

imagesc(QQ);


%%

for xx=1:(N_End-N_Start+1)
    if xx<9
        Background=Image_Stack(:,:,11);
    else
        Background=Image_Stack(:,:,1);
    end
    Image_Stack_Sub(:,:,xx)=Image_Stack(:,:,xx)-Background;
    disp(xx);
end

%%
NN=5;

imagesc(Image_Stack_Sub(:,:,NN));
colormap(gray);
caxis([0 500]);
axis equal
%%
for N=1:1

N=7;
Image=Image_Stack_Sub(:,:,N);
%Image=double(Image(ROI(1):ROI(2),ROI(3):ROI(4)));
imagesc(Image);
colormap(gray);
axis equal
%
Y=repmat([1:size(Image,1)]',1,size(Image,2));
X=repmat(1:size(Image,2),size(Image,1),1);


% Transfer to Polar Cord.
Wavelength=0.830;   %micron
X_Center=960;
Y_Center=790;
R_consider=1000;
NNN=2000;
LF=12;
HF=335;

R=repmat([R_consider/NNN:R_consider/NNN:R_consider],NNN,1);
THETA=repmat([2*pi/NNN:2*pi/NNN:2*pi]',1,NNN);

X_wish=R.*cos(THETA)+X_Center;
Y_wish=R.*sin(THETA)+Y_Center;

Image_Polar_Original=interp2(X,Y,Image,X_wish,Y_wish);
Image_Polar_Original(isnan(Image_Polar_Original))=0;
imagesc(Image_Polar_Original);

%
Image_Polar_RFFT=fft(Image_Polar_Original,[],2);
imagesc(real(Image_Polar_RFFT));


plot(real(Image_Polar_RFFT(1000,:)));

Image_Polar_RFFT(:,1:LF)=0;
Image_Polar_RFFT(:,HF:end)=0;
%Image_Polar_RFFT(:,(size(Image_Polar_RFFT,2)/2+1):end)=0;
Image_Polar_Filtered=ifft(Image_Polar_RFFT,[],2);
%Image_Polar_Filtered_Envelope=abs(ifft(Image_Polar_RFFT,[],2));
Image_Polar_Filtered_Envelope=abs(ifft(Image_Polar_RFFT,[],2));
imagesc(real(Image_Polar_Filtered));
%

%
WindowSize=30;


filter=fspecial('gaussian',WindowSize,round(WindowSize/2));
for p=1:size(filter,1)
    if p~=round(size(filter,1)/2)
        filter(:,p)=repmat(0,size(filter,2),1);
    end
end
%imagesc(filter);

filter=filter/sum(sum(filter));
%Image_Polar=conv2(Image_Polar_Original,filter,'same');
Image_Polar=conv2(real(Image_Polar_Filtered),filter,'same');

Image_Polar=max(max(Image_Polar))-Image_Polar;



imagesc(Image_Polar);
axis equal

% select特定等高線
Rough_R=1080;
first_tolerance=20;
number_of_line=20;
tolerance=1;
Next_Start_Index=1;
show_line=1;
clear Max_line_index
for q=1:number_of_line
    for p=1:NNN
        if (q==1)&&(p==1)
            [max_value max_index]=max(Image_Polar(p,(Rough_R-first_tolerance):(Rough_R+first_tolerance)));
            Max_line_index(p,q)=max_index+Rough_R-first_tolerance-1;
            First_index_record=Max_line_index(p,q);
        elseif (p==1)
            [max_value max_index]=max(Image_Polar(p,(Next_Start_Index-tolerance):(Next_Start_Index+tolerance)));
            Max_line_index(p,q)=max_index+Next_Start_Index-tolerance-1;
            
        else
            [max_value max_index]=max(Image_Polar(p,(Max_line_index(p-1,q)-tolerance):(Max_line_index(p-1,q)+tolerance)));
            Max_line_index(p,q)=max_index+Max_line_index(p-1,q)-tolerance-1;

        end
    end
    [next_min_value next_min_index]=findpeaks(max(max(Image_Polar))-Image_Polar(1,(Max_line_index(1,q)):end),'NPEAKS',1);
    [next_man_value next_max_index]=findpeaks(Image_Polar(1,(Max_line_index(1,q)+next_min_index):end),'NPEAKS',1);
    Next_Start_Index=next_min_index+next_max_index+Max_line_index(1,q)-1;
end
Last_index_record=Max_line_index(p,q);

    %findpeaks(Image_Polar(1,Max_line_index(1,q):end),'NPEAKS',1)+Max_line_index(1,q)-1;
if show_line==1
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
    %xlim([First_index_record Last_index_record]);
end
%%
    % back to the original RTheta base
basic_relative_height_array=-1*Wavelength/2*(1:number_of_line);
for p=1:NNN
    Z_RTheta(p,1:NNN)=interp1(Max_line_index(p,:),basic_relative_height_array,1:NNN);
end
%Z_RTheta(isnan(Z_RTheta))=0;


R_XY=((X-X_Center).^2+(Y-Y_Center).^2).^0.5;
THETA_XY=atan2((Y-Y_Center),(X-X_Center))+pi;

Z_XY(:,:,xx)=interp2(R,THETA,Z_RTheta,R_XY,THETA_XY);
imagesc(Z_XY(:,:,xx));
axis equal
colormap(gray);

%clear Max_line_index
end

Z_XY_backup=Z_XY;

TFmap_1=not(isnan(Z_XY(:,:,1)));
TFmap_2=not(isnan(Z_XY(:,:,2)));
TFmap_3=not(isnan(Z_XY(:,:,3)));
TFmap_4=not(isnan(Z_XY(:,:,4)));
TFmap_5=not(isnan(Z_XY(:,:,5)));
TFmap_6=not(isnan(Z_XY(:,:,6)));
TFmap_7=not(isnan(Z_XY(:,:,7)));
TFmap_8=not(isnan(Z_XY(:,:,8)));
TFmap_9=not(isnan(Z_XY(:,:,9)));
%%
NN=9;
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

for N=2:9
    if N==2
        Sample_index=[960 969];
    elseif N==3
        Sample_index=[1094 969];
    elseif N==4
        Sample_index=[1130 969];
    elseif N==5
        Sample_index=[1240 969];
    elseif N==6
        Sample_index=[1305 969];
    elseif N==7
        Sample_index=[1376 969];
    elseif N==8
        Sample_index=[1462 969];
    elseif N==9
        Sample_index=[1533 969];
    end
    Z_XY(:,:,N)=Z_XY(Sample_index(1),Sample_index(2),N-1)-Z_XY(Sample_index(1),Sample_index(2),N)+Z_XY(:,:,N);

end

plot((1:2040)*Lateral_Spacing,Z_XY(1:2040,969,1),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,2),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,3),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,4),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,5),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,6),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,7),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,8),(1:2040)*Lateral_Spacing,Z_XY(1:2040,969,9),'LineWidth',2);
xlabel('(micron)');
ylabel('(micron)');
xlim([0 5200]);
%%
Z_XY(isnan(Z_XY))=0;
Z_allsumup=sum(Z_XY,3);
imagesc(Z_allsumup);

TFmap_total=TFmap_1+TFmap_2+TFmap_3+TFmap_4+TFmap_5+TFmap_6+TFmap_7+TFmap_8+TFmap_9;
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

Pixel_size=5.5;
M=0.6;
Lateral_Spacing=Pixel_size*M;
k2=0.03;


R_Reference=10.3;%9;%7.79;          %mm
conic_Ref=0;

R_original=Lateral_Spacing*(((X-X(1,X_Center)).^2+(Y-Y(Y_Center,1)).^2).^0.5);
%imagesc(R_original);


R=R_original./(1+k2*0.000001*R_original.^2);
 plot(R(:,960),R_original(:,960));

ylabel('r'' (micron)');
xlabel('r (micron)');

%
Z_ref=((R.^2)/(R_Reference*1000))./(1+(1-(1+conic_Ref).*(R/(R_Reference*1000)).^2).^0.5);

Z_ref(Z_mod==0)=0;

% Axial Curvature Calculation
%
Z_offset=-2;
Z_total=Z_ref+Z_mod+Z_offset;
%Sin_Theta=R./(R.^2+Z.^2);
%Theta=asin(Sin_Theta);
%plot(Z_total(:,960));
%ylim([0 150]);
%xlim([500 1000]);
%
Theta=atan(R./abs(Z_total));
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

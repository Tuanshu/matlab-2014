clear all

cd('D:\Users\TuanShu\140127\lens_2\');

Stage_speed=30;  %micron/sec
Sampling_rate=350;   %Hz

Frame_Axial_Spacing=Stage_speed/Sampling_rate*5/3;  %micron

Objective_Focal_Length=4/0.85;   %cm
Porjection_Lens_Focal_Length=50;    %cm
Pixel_Size=5;       %micron

%Lateral_Spacing=Pixel_Size/Porjection_Lens_Focal_Length*Objective_Focal_Length;
Lateral_Spacing=Pixel_Size/20*5;

image_index=1:1990;
ROI=[1 508;1 271];      %up, down, left, right
Image_Stack(length(image_index),1:(ROI(1,2)-ROI(1,1)+1),1:(ROI(2,2)-ROI(2,1)+1))=0;
Low_Pass=600;   %pixel
High_Pass=800;
for p=1:length(image_index)
    Image=dlmread(sprintf('140127_lens_2_%i',image_index(p)));    
    %Image=imread(sprintf('1.jpg',image_index(p)));    
    Image_Stack(p,:,:)=Image(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2),1);
    disp(p);
end

FFT_Stack=fft(Image_Stack,[],1);
FFT_Stack(1:Low_Pass,:,:)=0;
FFT_Stack(round(length(FFT_Stack)/2):end,:,:)=0;
Image_Stack_New=real(ifft(FFT_Stack,[],1));
%plot(real(FFT_Stack(:,200,150)));

[max_value max_index]=max(Image_Stack_New,[],1);

MaxValue(:,:)=max_value(1,:,:);

%plot(real(FFT_Stack(:,300,400)));

Height(:,:)=Frame_Axial_Spacing*max_index(1,:,:);
Height=max(max(Height))-Height;
[maxvalue1D maxindex1D]=max(Height(:));
Xgrid(1:size(Height,1),1:size(Height,2))=0;
Ygrid(1:size(Height,1),1:size(Height,2))=0;

for p=1:size(Height,2)
    Xgrid(:,p)=(p-1)*Lateral_Spacing;
end
for q=1:size(Height,1)
    Ygrid(q,:)=(q-1)*Lateral_Spacing;
end

Xmax=Xgrid(maxindex1D);
Ymax=Ygrid(maxindex1D);
FitSurface= fittype( @(c,r,a, b, x, y) c+(r^2-(x-a).^2-(y-b).^2).^0.5, 'independent', {'x', 'y'},'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2
                                                                                                                        %z=C+(R^2-(x-A)^.2-(y-B)^.2).^0.5
                                                                                                                        
Starting_X=1450;                                                                                                                                                                                                                              
Starting_Y=1250;
R_Weight=450;
Weight_Func=MaxValue/max(max(MaxValue));
Weight_Func(:,:)=1;
Weight_Func((Xgrid-Starting_X).^2+(Ygrid-Starting_Y).^2>R_Weight^2)=0;
FitPara=fit([Xgrid(:),Ygrid(:)],Height(:),FitSurface,'Weight',Weight_Func(:),'StartPoint',[-15000,15000,Starting_X,Starting_Y]);

CheckCurve=FitPara.c+((FitPara.r)^2-(Xgrid-FitPara.a).^2-(Ygrid-FitPara.b).^2).^0.5;

imagesc(MaxValue/max(max(MaxValue)),'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal


imagesc(Weight_Func,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(MaxValue,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal

imagesc(CheckCurve,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,1)-1));
xlabel('(Micron)');
ylabel('(Micron)');
axis equal

plot((1:(size(Image_Stack_New,1)))*Frame_Axial_Spacing,real(Image_Stack_New(:,240,125)));
xlabel('Axial Position (Micron)');
ylabel('Amplitude (a.u.)');

Height2=Height;
Height2((Xgrid-Starting_X).^2+(Ygrid-Starting_Y).^2>R_Weight^2)=0;
plot(0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),Height2(250,:),0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),CheckCurve(250,:));
xlabel('(Micron)');
ylabel('(Micron)');
legend('Measured height','Fitting curve')


imagesc(Height);%,'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Height,1)-1));
xlabel('(Pixel)');
ylabel('(Pixel)');
axis equal
xlim([0 271]);
caxis([70 115]);
%%
Fitted_Curvature=FitPara.r %(micron)


%%
Binning_Factor_Height=32;
Pixel_Size_New=Pixel_Size*Binning_Factor_Height;

Height_Temp=Height;
[m,n]=size(Height_Temp);

Height_Temp=sum(reshape(Height_Temp,Binning_Factor_Height,[]),1);
Height_Temp=reshape(Height_Temp,m/Binning_Factor_Height,[]).';

Height_Temp=sum( reshape(Height_Temp,Binning_Factor_Height,[]) ,1);
Height_Temp=reshape(Height_Temp,n/Binning_Factor_Height,[]).';


Height_Temp=Height_Temp/Binning_Factor_Height^2;

Height_Binned=Height_Temp;

imagesc(Height_Temp);
%%

dz_dx=diff(Height_Binned,1,1)/Pixel_Size_New;
dz_dx(size(dz_dx,1)+1,:)=dz_dx(size(dz_dx,1),:);
d2z_dx2=diff(Height_Binned,2,1)/Pixel_Size_New^2;
d2z_dx2(size(d2z_dx2,1)+1,:)=d2z_dx2(size(d2z_dx2,1),:);
d2z_dx2(size(d2z_dx2,1)+1,:)=d2z_dx2(size(d2z_dx2,1),:);

dz_dy=diff(Height_Binned,1,2)/Pixel_Size_New;
dz_dy(:,size(dz_dy,2)+1)=dz_dy(:,size(dz_dy,2));

d2z_dy2=diff(Height_Binned,2,2)/Pixel_Size_New^2;
d2z_dy2(:,size(d2z_dy2,2)+1)=d2z_dy2(:,size(d2z_dy2,2));
d2z_dy2(:,size(d2z_dy2,2)+1)=d2z_dy2(:,size(d2z_dy2,2));

%%
Kx=d2z_dx2./(1+(dz_dx).^2).^(3/2);
Ky=d2z_dy2./(1+(dz_dy).^2).^(3/2);

%%

Rx=1./Kx/1000;  %mm
Ry=1./Ky/1000;  %mm
imagesc(Rx);
caxis([0 20]);

%%

N=143;

imagetemp(:,:)=Image_Stack(N,:,:);
imagesc(imagetemp);
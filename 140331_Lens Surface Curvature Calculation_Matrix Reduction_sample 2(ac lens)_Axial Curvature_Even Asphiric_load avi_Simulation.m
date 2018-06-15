clear all

d_Position=0.56/2/49;
Lateral_Spacing=122/62;

%% Height generation
cd('D:\Users\TuanShu\140304\lens_3\');
Averaging_Factor=1;
Binning_Factor=1;

SPF=25;
LPF=55;   %pixel
%% Height generation
cd('E:\Users\TuanShu\140331_Crystalvue');
%MOVIE=aviread('140331_30x-30x.avi');
%%
%imagesc(MOVIE(191).cdata);
%imagetest=MOVIE(191).cdata;
%image_index=1:1840;

ROI=[1 344;1 768];      %up, down, left, right


X(1:ROI(1,2),1:ROI(2,2))=0;
Y(1:ROI(1,2),1:ROI(2,2))=0;

for p=1:size(X,2)
    X(:,p)=p.*Lateral_Spacing;
end
for q=1:size(Y,1)
    Y(q,:)=q.*Lateral_Spacing;
end




%Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask)*(-1);
%Z=((Surface_profile-max(max(Surface_profile))).*d_Position.*Mask_new)*(-1);
%Z(isnan(Z))=0;




%% Axial Curvature Calculation



X_offset=430;
Y_offset=180;

R=((X-X(1,X_offset)).^2+(Y-Y(Y_offset,1)).^2).^0.5;
imagesc(R);
%%
R_mm=R/1000;

K_1=-87.22806;
C_inverse_1=-14.57544;
C_1=1/C_inverse_1;
coef2_1=0;
coef4_1=-9.15744E-4;
coef6_1=8.643373E-5;
coef8_1=0;
coef10_1=0;

R_Z2=-14.66;

Z1=((C_1.*(R_mm.^2))./(1+(1-(1+K_1).*(C_1.^2).*(R_mm.^2)).^0.5)+coef2_1.*(R_mm.^2)+coef4_1.*(R_mm.^4)+coef6_1.*(R_mm.^6)+coef8_1.*(R_mm.^8)+coef10_1.*(R_mm.^10));
imagesc(Z1);

K_2=0;
C_inverse_2=-19.9233;
C_2=1/C_inverse_2;
coef2_2=0;
coef4_2=0.0042371;
coef6_2=-0.00063484;
coef8_2=3.46526E-5;
coef10_2=0;

Z2=((C_2.*(R_mm.^2))./(1+(1-(1+K_2).*(C_2.^2).*(R_mm.^2)).^0.5)+coef2_2.*(R_mm.^2)+coef4_2.*(R_mm.^4)+coef6_2.*(R_mm.^6)+coef8_2.*(R_mm.^8)+coef10_2.*(R_mm.^10));
imagesc(Z2);




%%
Z=(Z2)*1000;
%Sin_Theta=R./(R.^2+Z.^2);
%Theta=asin(Sin_Theta);

Theta=atan(R./abs(Z));
Sin_pi_minus_2Theta=sin(pi-2*Theta);
Axial_Radius_of_Curvature=(R./Sin_pi_minus_2Theta)/1000;
Diopter=(1.3375-1)*(1000)./Axial_Radius_of_Curvature;

imagesc(R,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
imagesc(Axial_Radius_of_Curvature,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
caxis([0 30]);
axis equal
xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);

colorbar
xlabel('(micron)');
ylabel('(micron)');


%imagesc(Diopter,'xdata',Lateral_Spacing:Lateral_Spacing:(size(Z,2)*Lateral_Spacing),'ydata',Lateral_Spacing:Lateral_Spacing:(size(Z,1)*Lateral_Spacing));
%axis equal
%xlim([Lateral_Spacing size(Z,2)*Lateral_Spacing]);
%ylim([Lateral_Spacing size(Z,1)*Lateral_Spacing]);
%caxis([0 50]);

%colorbar
%xlabel('(micron)');
%ylabel('(micron)');
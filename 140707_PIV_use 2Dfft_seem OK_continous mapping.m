clear all

%% Note: ¶Âªº¤~¬Ocell

interrogation_size=20;

Total_Size_X=480;
Total_Size_Y=640;

Pixel_Spacing=0.45;     %micron

Frame_Rate=250;%*526/2000;         %fps
Normal_Velocity=Pixel_Spacing*Frame_Rate;  %micron/second

cd('D:\Users\TuanShu\140626_RBC\OPSI\tif\');
File_list=dir;
%MOVIE=aviread('red blood cell in tissue.avi');

%Test=MOVIE(205).cdata(:,:,1)-MOVIE(205).cdata(:,:,3);

%%
Threshold=20;

Blur_Window_Size=10;


filter=fspecial('gaussian',Blur_Window_Size,round(Blur_Window_Size/2));
filter=filter/sum(sum(filter));

for p=1:length(File_list)-2
    Image(:,:,p)=conv2(double(imread(File_list(p+2).name)),filter,'same');
    disp(p);
end
%%
%LP=3;
%FFT_Image=fft(Image,[],3);

%FFT_Image(:,:,1:LP)=0;
%FFT_Image(:,:,(size(FFT_Image,3)/2+1):size(FFT_Image,3))=0;

%Image_New=abs(ifft(FFT_Image,[],3));

%%
Background=mean(Image,3);
Image=repmat(Background,[1 1 size(Image,3)])-Image;

Image(Image<Threshold)=0;
%%


N=38;
%Image_New(Image_New<Threshold)=0;

imagesc(Image(:,:,N));

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
%% Second Blur
Blur_Window_Size_2=10;

filter_2=fspecial('gaussian',Blur_Window_Size_2,round(Blur_Window_Size_2/2));
filter_2=filter_2/sum(sum(filter_2));


for p=1:size(Image,3)
    Image(:,:,p)=conv2(Image(:,:,p),filter_2,'same');
    disp(p);
end
Image(Image<Threshold)=0;

%%
%
N=49;
%Image_New(Image_New<Threshold)=0;

imagesc(Image(:,:,N));

%% Output avi


%% All

clear Current_Interrogation_1 Current_Interrogation_2 Current_Interrogation_1_ZP Current_Interrogation_2_ZP V_CM_X V_CM_Y Sum_V_CM_X Sum_V_CM_Y
original_interrogation_size=20;
interpolation=5;
interrogation_size=interpolation*original_interrogation_size;

X_grid_old=repmat((1:original_interrogation_size)'*interpolation,1,original_interrogation_size);
Y_grid_old=repmat((1:original_interrogation_size)*interpolation,original_interrogation_size,1);

X_grid_new=repmat((1:interrogation_size)',1,interrogation_size);
Y_grid_new=repmat((1:interrogation_size),interrogation_size,1);

Start_Frame=1;
Number_of_Frame_Considered=20;
Frame_Resolution=25;


Offset_X_Start=160;
Offset_Y_Start=80;
Resolution=2;

N_X=100;
N_Y=100;

%Offset_X_2=Offset_X+0;
%Offset_Y_2=Offset_Y+0;


Zero_Padding_N=1;           %N each side

Current_Interrogation_1_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;
Current_Interrogation_2_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;

cordinate_X=repmat((1:N_X)',1,N_Y);
cordinate_Y=repmat((1:N_Y),N_X,1);

V_CM_X(1:N_X,1:N_Y)=0;
Sum_V_CM_X(1:N_X,1:N_Y)=0;
V_CM_Y(1:N_X,1:N_Y)=0;
Sum_V_CM_Y(1:N_X,1:N_Y)=0;
for w=1:Number_of_Frame_Considered
    
    Frame_N_1=Start_Frame+Frame_Resolution*(w-1);%124;
    Frame_N_2=Frame_N_1+1;%124;

    for p=1:N_X
        for q=1:N_Y
            Current_Interrogation_1=interp2(Y_grid_old,X_grid_old,Image((Offset_X_Start+(p-1)*Resolution):(original_interrogation_size+Offset_X_Start+(p-1)*Resolution-1),(Offset_Y_Start+(q-1)*Resolution):(original_interrogation_size+Offset_Y_Start+(q-1)*Resolution-1),Frame_N_1),Y_grid_new,X_grid_new);
            Current_Interrogation_1(isnan(Current_Interrogation_1))=0;
            Current_Interrogation_2=interp2(Y_grid_old,X_grid_old,Image((Offset_X_Start+(p-1)*Resolution):(original_interrogation_size+Offset_X_Start+(p-1)*Resolution-1),(Offset_Y_Start+(q-1)*Resolution):(original_interrogation_size+Offset_Y_Start+(q-1)*Resolution-1),Frame_N_2),Y_grid_new,X_grid_new);
            Current_Interrogation_2(isnan(Current_Interrogation_2))=0;
            Current_Interrogation_1_ZP((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size))=Current_Interrogation_1;
            Current_Interrogation_2_ZP((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size))=Current_Interrogation_2;


            FFT_1=fft2(Current_Interrogation_1_ZP);
            FFT_2=fft2(Current_Interrogation_2_ZP);
            V_Temp=ifftshift(ifftshift(ifft2(conj(FFT_1).*FFT_2),1),2);

            V=V_Temp((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size));

            %subplot(2,2,1)
            %imagesc(Current_Interrogation_1);
            %subplot(2,2,2)
            %imagesc(Current_Interrogation_2);
            %subplot(2,2,3)
            %imagesc((V_Temp));
            %subplot(2,2,4)
            %imagesc((V_Temp));


            %V_X_CM_Y=sum(sum(V_X.*Y_grid))/sum(sum(V_X));
            %V_Y_CM_X=sum(sum(V_Y.*X_grid))/sum(sum(V_Y));
            V_CM_X(p,q)=sum(sum(V.*X_grid_new))/sum(sum(V))-1-interrogation_size/2;
            V_CM_Y(p,q)=sum(sum(V.*Y_grid_new))/sum(sum(V))-1-interrogation_size/2;


            %fprintf('(V_CM_X,V_CM_Y)=(%f,%f)\n',V_CM_X(p,q),V_CM_Y(p,q));

        end
        %disp(p);

    end
    V_CM_X(isnan(V_CM_X))=0;
    V_CM_Y(isnan(V_CM_Y))=0;
    Sum_V_CM_X=Sum_V_CM_X+V_CM_X;
    Sum_V_CM_Y=Sum_V_CM_Y+V_CM_Y;
    disp(w);
end

Ave_V_CM_X=Sum_V_CM_X/Number_of_Frame_Considered;
Ave_V_CM_Y=Sum_V_CM_Y/Number_of_Frame_Considered;
%h=quiver3(X_resam,Y_resam,Z_resam,Xr_XY_resam,Yr_XY_resam,Zr_XY_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
%adjust_quiver_arrowhead_size(h,2);
%
%Current_Image(:,:)=Image(:,:,320);
%imagesc(Current_Image);

%axis equal
%xlim([1 648]);
%ylim([1 488]);
%        set(gca, 'XTick', []);
%        set(gca, 'YTick', []);
%        set(gca,'XColor','white');
%        set(gca,'YColor','white');
%hold on

%Ave_V_CM_Y(8,26)=0;
%Ave_V_CM_X(8,26)=0;
%Ave_V_CM_Y(13,12)=0;
%Ave_V_CM_X(13,12)=0;
%h=quiver(cordinate_Y+original_interrogation_size/2,cordinate_X+original_interrogation_size/2,Ave_V_CM_Y,Ave_V_CM_X,'linewidth',2,'color',[1 1 1]);%'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
%axis equal
%xlim([1 648]);
%ylim([1 488]);

%hold off

%

Velocity=((Ave_V_CM_Y.^2+Ave_V_CM_X.^2).^0.5)*Normal_Velocity;

imagesc(Velocity);

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        colorbar
Velocity_N=Velocity;
%Velocity_N(8,26)=0;
%Velocity_N(13,12)=0;

imagesc(Velocity_N/2,'xdata',(Offset_X_Start+(1:N_X)*Resolution),'ydata',(Offset_Y_Start+(1:N_Y)*Resolution));

        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca,'XColor','white');
        set(gca,'YColor','white');
        colorbar
        axis equal 
        
%xlim([1 32]);
%ylim([1 24]);
caxis([0 500]);


%% backup
V_map=Velocity_N/2;
imagesc(V_map);

Blur_Window_Size_3=3;

filter_3=fspecial('gaussian',Blur_Window_Size_3,round(Blur_Window_Size_3/2));
filter_3=filter_3/sum(sum(filter_3));

V_map(V_map>500)=200;

V_map=conv2(V_map,filter_3,'same');

imagesc(V_map,'ydata',(Offset_X_Start+(1:N_X)*Resolution),'xdata',(Offset_Y_Start+(1:N_Y)*Resolution));
xlabel('(Pixel)');
ylabel('(Pixel)');
colormap(jet);
colorbar
caxis([0 500]);
%%
imagesc(double(imread(File_list(p+2).name)),'xdata',1:648,'ydata',1:488);
colormap(gray);

hold on



imagesc(V_map,'ydata',(Offset_X_Start+(1:N_X)*Resolution),'xdata',(Offset_Y_Start+(1:N_Y)*Resolution));
colormap(jet);
caxis([0 500]);

hold off


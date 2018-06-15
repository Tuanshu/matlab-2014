clear all

%% Note: 黑的才是cell

interrogation_size=20;

Total_Size_X=480;
Total_Size_Y=640;



cd('D:\Users\TuanShu\140626_RBC\OPSI\tif\');
File_list=dir;
%MOVIE=aviread('red blood cell in tissue.avi');

%Test=MOVIE(205).cdata(:,:,1)-MOVIE(205).cdata(:,:,3);

%%
Threshold=70;
for p=1:length(File_list)-2
    Image(:,:,p)=double(imread(File_list(p+2).name));
    disp(p);
end
%
%N=282;
%colormap(gray);
%imagesc(Image(:,:,N));

%%
Background=mean(Image,3);
Image=repmat(Background,[1 1 size(Image,3)])-Image;

Image(Image<Threshold)=0;


%%
N=2;

imagesc(Image(:,:,N));

%%
%% Single
clear Current_Interrogation_1 Current_Interrogation_2 Current_Interrogation_1_ZP Current_Interrogation_2_ZP V_X_CM_X V_Y_CM_Y
original_interrogation_size=20;
interpolation=2;
interrogation_size=interpolation*original_interrogation_size;

X_grid_old=repmat((1:original_interrogation_size)'*interpolation,1,original_interrogation_size);
Y_grid_old=repmat((1:original_interrogation_size)*interpolation,original_interrogation_size,1);


X_grid_new=repmat((1:interrogation_size)',1,interrogation_size);
Y_grid_new=repmat((1:interrogation_size),interrogation_size,1);

p=10;
q=10;

Offset_X=18;
Offset_Y=-13;

Offset_X_2=Offset_X+5;
Offset_Y_2=Offset_Y+0;

Frame_N_1=3;%124;
Frame_N_2=Frame_N_1+0;%124;

Zero_Padding_N=1;           %N each side

Current_Interrogation_1_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;
Current_Interrogation_2_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;


cordinate_X=repmat((1:(Total_Size_X/interrogation_size))',1,(Total_Size_Y/interrogation_size));
cordinate_Y=repmat((1:(Total_Size_Y/interrogation_size)),(Total_Size_X/interrogation_size),1);

V_X_CM_X(1:(Total_Size_X/interrogation_size),1:(Total_Size_Y/interrogation_size))=0;
V_Y_CM_Y(1:(Total_Size_X/interrogation_size),1:(Total_Size_Y/interrogation_size))=0;

        Current_Interrogation_1=interp2(Y_grid_old,X_grid_old,Image((1+(p-1)*original_interrogation_size+Offset_X):((p)*original_interrogation_size+Offset_X),(1+(q-1)*original_interrogation_size+Offset_Y):((q)*original_interrogation_size+Offset_Y),Frame_N_1),Y_grid_new,X_grid_new);
        Current_Interrogation_1(isnan(Current_Interrogation_1))=0;
        Current_Interrogation_2=interp2(Y_grid_old,X_grid_old,Image((1+(p-1)*original_interrogation_size+Offset_X_2):((p)*original_interrogation_size+Offset_X_2),(1+(q-1)*original_interrogation_size+Offset_Y_2):((q)*original_interrogation_size+Offset_Y_2),Frame_N_2),Y_grid_new,X_grid_new);
        Current_Interrogation_2(isnan(Current_Interrogation_2))=0;
        Current_Interrogation_1_ZP((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size))=Current_Interrogation_1;
        Current_Interrogation_2_ZP((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size))=Current_Interrogation_2;


        
        FFT_1_X=fft(Current_Interrogation_1_ZP,[],1);
        FFT_2_X=fft(Current_Interrogation_2_ZP,[],1);
        FFT_1_Y=fft(Current_Interrogation_1_ZP,[],2);
        FFT_2_Y=fft(Current_Interrogation_2_ZP,[],2);
        
        V_X_Temp=ifftshift(ifft(conj(FFT_1_X).*FFT_2_X,[],1),1);
        V_Y_Temp=ifftshift(ifft(conj(FFT_1_Y).*FFT_2_Y,[],2),2);
        
        V_X=V_X_Temp((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size));
        V_Y=V_Y_Temp((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size));

        subplot(2,2,1)
        imagesc(Current_Interrogation_1);
        subplot(2,2,2)
        imagesc(Current_Interrogation_2);
        subplot(2,2,3)
        imagesc(V_X);
        subplot(2,2,4)
        imagesc(V_Y);
        
        
        V_X_CM_X(p,q)=sum(sum(V_X.*X_grid_new))/sum(sum(V_X))-1-interrogation_size/2;
        %V_X_CM_Y=sum(sum(V_X.*Y_grid))/sum(sum(V_X));
        %V_Y_CM_X=sum(sum(V_Y.*X_grid))/sum(sum(V_Y));
        V_Y_CM_Y(p,q)=sum(sum(V_Y.*Y_grid_new))/sum(sum(V_Y))-1-interrogation_size/2;
        
        fprintf('(V_X_CM_X,V_Y_CM_Y)=(%f,%f)\n',V_X_CM_X(p,q),V_Y_CM_Y(p,q));
        disp(p);
      %%  
        % to calculate the center of mass of V_X and V_Y (Q:對V_X V_Y而言,
        % 另外一個dim的center of mass的移動意味著
        X_grid=repmat((1:interrogation_size)',1,interrogation_size)-ceil(interrogation_size/2)-1;
        Y_grid=repmat((1:interrogation_size),interrogation_size,1)-ceil(interrogation_size/2)-1;
        
        V_X_CM_X(p,q)=sum(sum(V_X.*X_grid))/sum(sum(V_X));
        %V_X_CM_Y=sum(sum(V_X.*Y_grid))/sum(sum(V_X));
        %V_Y_CM_X=sum(sum(V_Y.*X_grid))/sum(sum(V_Y));
        V_Y_CM_Y(p,q)=sum(sum(V_Y.*Y_grid))/sum(sum(V_Y));
        
        %fprintf('(V_X_CM_X,V_Y_CM_Y)=(%f,%f)\n',V_X_CM_X,V_Y_CM_Y);
        %fprintf('(V_Y_CM_X,V_X_CM_Y)=(%f,%f)\n',V_Y_CM_X,V_X_CM_Y);
        disp(p);

%h=quiver3(X_resam,Y_resam,Z_resam,Xr_XY_resam,Yr_XY_resam,Zr_XY_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
%adjust_quiver_arrowhead_size(h,2);
%
Current_Image(:,:)=Image(:,:,Frame_N_1);
imagesc(Current_Image);

xlabel('(Micron)');
ylabel('(Micron)');
zlabel('(Micron)');

hold on
h=quiver(cordinate_Y*interrogation_size,cordinate_X*interrogation_size,V_Y_CM_Y,V_X_CM_X,'linewidth',2,'color',[1 1 1]);%'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
hold off




%% All
clear Current_Interrogation_1 Current_Interrogation_2 Current_Interrogation_1_ZP Current_Interrogation_2_ZP V_X_CM_X V_Y_CM_Y
interrogation_size=10;
Offset_X=0;
Offset_Y=0;

Offset_X_2=Offset_X+0;
Offset_Y_2=Offset_Y+0;

Frame_N_1=1;%124;
Frame_N_2=Frame_N_1+1;%124;

Zero_Padding_N=1;           %N each side

Current_Interrogation_1_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;
Current_Interrogation_2_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;


cordinate_X=repmat((1:(Total_Size_X/interrogation_size))',1,(Total_Size_Y/interrogation_size));
cordinate_Y=repmat((1:(Total_Size_Y/interrogation_size)),(Total_Size_X/interrogation_size),1);

V_X_CM_X(1:(Total_Size_X/interrogation_size),1:(Total_Size_Y/interrogation_size))=0;
V_Y_CM_Y(1:(Total_Size_X/interrogation_size),1:(Total_Size_Y/interrogation_size))=0;
for p=1:(Total_Size_X/interrogation_size)
    for q=1:(Total_Size_Y/interrogation_size)
        Current_Interrogation_1=Image((1+(p-1)*interrogation_size+Offset_X):((p)*interrogation_size+Offset_X),(1+(q-1)*interrogation_size+Offset_Y):((q)*interrogation_size+Offset_Y),Frame_N_1);
        Current_Interrogation_2=Image((1+(p-1)*interrogation_size+Offset_X_2):((p)*interrogation_size+Offset_X_2),(1+(q-1)*interrogation_size+Offset_Y_2):((q)*interrogation_size+Offset_Y_2),Frame_N_2);

        Current_Interrogation_1_ZP((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size))=Image((1+(p-1)*interrogation_size+Offset_X):((p)*interrogation_size+Offset_X),(1+(q-1)*interrogation_size+Offset_Y):((q)*interrogation_size+Offset_Y),Frame_N_1);
        Current_Interrogation_2_ZP((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size))=Image((1+(p-1)*interrogation_size+Offset_X_2):((p)*interrogation_size+Offset_X_2),(1+(q-1)*interrogation_size+Offset_Y_2):((q)*interrogation_size+Offset_Y_2),Frame_N_2);


        
        FFT_1_X=fft(Current_Interrogation_1_ZP,[],1);
        FFT_2_X=fft(Current_Interrogation_2_ZP,[],1);
        FFT_1_Y=fft(Current_Interrogation_1_ZP,[],2);
        FFT_2_Y=fft(Current_Interrogation_2_ZP,[],2);
        
        V_X_Temp=ifftshift(ifft(conj(FFT_1_X).*FFT_2_X,[],1),1);
        V_Y_Temp=ifftshift(ifft(conj(FFT_1_Y).*FFT_2_Y,[],2),2);
        
        V_X=V_X_Temp((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size));
        V_Y=V_Y_Temp((Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size),(Zero_Padding_N*interrogation_size+1):((Zero_Padding_N+1)*interrogation_size));

        %subplot(2,2,1)
        %imagesc(Current_Interrogation_1);
        %subplot(2,2,2)
        %imagesc(Current_Interrogation_2);
        %subplot(2,2,3)
        %imagesc(V_X);
        %subplot(2,2,4)
        %imagesc(V_Y);
        
        % to calculate the center of mass of V_X and V_Y (Q:對V_X V_Y而言,
        % 另外一個dim的center of mass的移動意味著
        X_grid=repmat((1:interrogation_size)',1,interrogation_size)-ceil(interrogation_size/2)-1;
        Y_grid=repmat((1:interrogation_size),interrogation_size,1)-ceil(interrogation_size/2)-1;
        
        V_X_CM_X(p,q)=sum(sum(V_X.*X_grid))/sum(sum(V_X));
        %V_X_CM_Y=sum(sum(V_X.*Y_grid))/sum(sum(V_X));
        %V_Y_CM_X=sum(sum(V_Y.*X_grid))/sum(sum(V_Y));
        V_Y_CM_Y(p,q)=sum(sum(V_Y.*Y_grid))/sum(sum(V_Y));
        
        %fprintf('(V_X_CM_X,V_Y_CM_Y)=(%f,%f)\n',V_X_CM_X,V_Y_CM_Y);
        %fprintf('(V_Y_CM_X,V_X_CM_Y)=(%f,%f)\n',V_Y_CM_X,V_X_CM_Y);
        disp(p);

    end
end

%h=quiver3(X_resam,Y_resam,Z_resam,Xr_XY_resam,Yr_XY_resam,Zr_XY_resam,'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
%adjust_quiver_arrowhead_size(h,2);
%
Current_Image(:,:)=Image(:,:,Frame_N_1);
imagesc(Current_Image);

xlabel('(Micron)');
ylabel('(Micron)');
zlabel('(Micron)');

hold on
h=quiver(cordinate_Y*interrogation_size,cordinate_X*interrogation_size,V_Y_CM_Y,V_X_CM_X,'linewidth',2,'color',[1 1 1]);%'linewidth',2,'color',[0 0 0],'xdata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,2)-1),'ydata',0:Lateral_Spacing:Lateral_Spacing*(size(Curvature_P1,1)-1));
hold off


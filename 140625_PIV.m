clear all

interrogation_size=20;

Total_Size_X=240;
Total_Size_Y=320;

Threshold=15;


cd('D:\Users\TuanShu\140625_PIV\');

MOVIE=aviread('red blood cell in tissue.avi');

Test=MOVIE(205).cdata(:,:,1)-MOVIE(205).cdata(:,:,3);

%
imagesc(Test);
%%
for p=1:800
    Image(:,:,p)=MOVIE(p).cdata(1:240,1:320,1)-MOVIE(p).cdata(1:240,1:320,3);
    disp(p);
end
Image(Image<Threshold)=0;

%% 
clear Current_Interrogation_1 Current_Interrogation_2 Current_Interrogation_1_ZP Current_Interrogation_2_ZP
interrogation_size=10;
Offset_X=16;
Offset_Y=25;

Offset_X_2=Offset_X+0;
Offset_Y_2=Offset_Y+0;

Frame_N_1=125;%124;
Frame_N_2=125;%124;

Zero_Padding_N=1;           %N each side

Current_Interrogation_1_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;
Current_Interrogation_2_ZP(1:(2*Zero_Padding_N+1)*interrogation_size,1:(2*Zero_Padding_N+1)*interrogation_size)=0;

p=14;
q=18;

%
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

        subplot(2,2,1)
        imagesc(Current_Interrogation_1);
        subplot(2,2,2)
        imagesc(Current_Interrogation_2);
        subplot(2,2,3)
        imagesc(V_X);
        subplot(2,2,4)
        imagesc(V_Y);
        
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

        end
end

%%
A(:,:)=Image(:,:,N);
B(:,:)=Image(:,:,N+10);

FFTA=fft(A,[],2);
FFTB=fft(B,[],2);

TEST=ifft(conj(FFTA).*FFTB,[],2);


subplot(3,1,1)

imagesc(A);

subplot(3,1,2)

imagesc(B);

subplot(3,1,3)

imagesc(TEST);

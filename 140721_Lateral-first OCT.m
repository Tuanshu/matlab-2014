clear all

Pixel_size=5.5;
M=0.6;
Lateral_Spacing=Pixel_size*M;
KK=350;
N_Start=360;
N_End=385;
% Height generation

cd('D:\Users\TuanShu\140721\12bit_front');
%% background

Image_Stack(1:2040,1:2048,1:(N_End-N_Start+1))=0;

for xx=1:(N_End-N_Start+1)

Image_Stack(:,:,xx)=double(dlmread(sprintf('140721_12bit_A_Eye_2_front_%d',N_Start+xx-1)));
disp(xx);
end


%%

NN=14;

QQ(:,:)=Image_Stack(:,:,NN);

imagesc(QQ);


%%

Background=mean(Image_Stack,3);

for xx=1:(N_End-N_Start-1)

Image_Stack_Sub(:,:,xx)=Image_Stack(:,:,xx+1)-Image_Stack(:,:,xx);
disp(xx);
end
%%
NN=15;

imagesc(Image_Stack_Sub(:,:,NN));
colormap(gray);
%caxis([-200 200]);
axis equal
%%
LP=50;

NN=7;
%%xx=NN;
for xx=1:size(Image_Stack_Sub,3)

FFT2(:,:,xx)=fft2(Image_Stack_Sub(:,:,xx));
disp(xx);

end

FFT2_Filtered=FFT2;
FFT2_Filtered(size(FFT2_Filtered,1)/2:end,:,:)=0;
FFT2_Filtered(1:LP,:,:)=0;

FFT2_Filtered(:,size(FFT2_Filtered,2)/2:end,:)=0;
FFT2_Filtered(:,1:LP,:)=0;


%
%NN=14;

%imagesc(fftshift(fftshift(real(FFT2_Filtered(:,:,NN)),2),1));
%caxis([-1E5 1E5]);

%
for xx=1:size(Image_Stack_Sub,3)
Image_New(:,:,xx)=ifft2(FFT2_Filtered(:,:,xx));
disp(xx);

end

%%
NN=16;
imagesc(abs(Image_New(:,:,NN)));

axis equal
%% output
for xx=1:size(Image_New,3)
    dlmwrite(sprintf('Output_%d',xx),Image_New(:,:,xx),'delimiter','\t','newline','pc','precision', '%.6f');
    disp(xx);
end


%%

V=0.1; %%mm/sec
Frame_rate=30;

Frame_Spacing=V/Frame_rate*1000;    %micton

[max_value max_index]=max(abs(Image_New),[],3);

Z=max_index*Frame_Spacing;
imagesc(Z);
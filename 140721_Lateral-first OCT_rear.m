clear all

Pixel_size=5.5;
M=0.6;
Lateral_Spacing=Pixel_size*M;
KK=350;
N_Start=540;
N_End=580;
% Height generation

cd('D:\Users\TuanShu\140721\12bit_rear');
%% background

Image_Stack(1:2040,1:2048,1:(N_End-N_Start+1))=0;

for xx=1:(N_End-N_Start+1)

Image_Stack(:,:,xx)=double(dlmread(sprintf('140721_12bit_A_Eye_3_rear_%d',N_Start+xx-1)));
disp(xx);
end


%%

NN=17;

QQ(:,:)=Image_Stack(:,:,NN);

imagesc(QQ);


%%

Background=mean(Image_Stack,3);

for xx=1:(N_End-N_Start-1)

Image_Stack_Sub(:,:,xx)=Image_Stack(:,:,xx+1)-Image_Stack(:,:,xx);
disp(xx);
end
%Image_Stack_Sub=Image_Stack-repmat(Background,[1 1 size(Image_Stack,3)]);
%%
NN=79;

imagesc(Image_Stack_Sub(:,:,NN));
colormap(gray);
%caxis([-200 200]);
axis equal
%%



for xx=1:size(Image_Stack_Sub,3)

FFT2(:,:,xx)=fft2(Image_Stack_Sub(:,:,xx));
disp(xx);

end

FFT2_Filtered=FFT2;
FFT2_Filtered(size(FFT2_Filtered,1)/2:end,:,:)=0;
FFT2_Filtered(1:2,:,:)=0;

FFT2_Filtered(:,size(FFT2_Filtered,2)/2:end,:)=0;
FFT2_Filtered(:,1:2,:)=0;


%%
NN=14;

imagesc(fftshift(fftshift(real(FFT2_Filtered(:,:,NN)),2),1));
caxis([-1E5 1E5]);

%%

Image_New=ifft2(FFT2_Filtered);
%%


NN=17;

imagesc(abs(Image_New(:,:,NN)));

axis equal

%%

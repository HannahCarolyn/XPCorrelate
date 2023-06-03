Map_PQ_before_KMC=imsegkmeans(uint16(Data_InputMap_before.RadonQuality),2,"NumAttempts",10);
Map_PQ_before_KMC=(im2double(Map_PQ_before_KMC)*255);

figure;
imagesc(Map_PQ_before_KMC)
axis image
colormap gray

PH2_before=Map_PH2_before;
figure;
I=imagesc(PH2_before);
colormap gray

IE=entropyfilt(PH2_before);
Eim=rescale(IE);
figure;
imshow(Eim)

BW1 = imbinarize(Eim,0.89);
figure;
imshow(BW1)

BW1double=double(BW1);

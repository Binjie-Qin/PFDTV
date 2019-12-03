image=imread('57811-Afbeelding3.jpg');
image=rgb2gray(image);
figure,imshow(image);
PC  = phasesymmono(image,4,3);
%[PC, or, ft, T]= phasecong2(image,4,3);
figure,imshow(PC);
imwrite(PC,'PS.jpg');
%save('PCMAP.txt','PC','-ASCII');
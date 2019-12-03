image=imread('origin.jpg');
%image=rgb2gray(image);
figure,imshow(image);
%PC  = phasesymmono(image,4,3);
[PC, or, ft, T]= phasecong2(image,4,3);
figure,imshow(PC);
save('PCMAP.txt','PC','-ASCII');
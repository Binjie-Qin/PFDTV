image=imread('origin.jpg');
%image=rgb2gray(image);
figure,imshow(image);
[rows,cols]=size(image);
IM = fft2(image);       
[radius, u1, u2] = FilterGrid(rows, cols);
radius(1,1) = 1;
sumPC=zeros(rows,cols);
[PCm, or, ft, Ts] = phasecongmono(image,5,5);
for s=1:15
     xishu=sqrt((pi*(4.^2.58)*(s.^4.16))/gamma(4.16));
     cauchy=xishu.*(radius.^1.58).*exp(-s*radius);
     H = (1i*u1 - u2)./radius; 
     IMF = IM.*cauchy; 
     f = real(ifft2(IMF)); 
     h = ifft2(IMF.*H);  
     
     h1 = real(h); 
     h2 = imag(h);
     An = sqrt(f.^2 + h1.^2 + h2.^2);
     odd=sqrt(h1.^2+h2.^2);
     even=sqrt(f.^2);
     
    PC=(abs(odd)-abs(even)-Ts)./(sqrt(f.^2 + h1.^2 + h2.^2)+0.0001);
    PC(PC<0)=0;
     
     sumPC=sumPC+PC;
end
PC=sumPC;
PC1 = mapminmax(PC, 0, 1);
figure,imshow(PC1);%imwrite(PC1,'25.jpg');
save('PCMAP.txt','PC','-ASCII');
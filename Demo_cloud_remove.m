clear,clc

load('L8_Jizzakh_R3_C10');
load('Mask');

[M,N,B1,B2] = size(data);

Ori_H=im2double(data);
Ori_H=reshape(Ori_H,M,N,B1*B2);
Mask=reshape(Mask,M,N,B1*B2);
Noi_H=im2double(Ori_H.*im2double(Mask));

[L2 E2 ]=HRPCA2D_MC(Noi_H,Mask,'gamma',[1,1,1,1],'maxIter',200,'rho',1.05,'lambda',2,'debug',1,'GT',Ori_H,'tol', 1e-6);

[psnr, ssim, fsim, ergas, msam] = MSIQA(L2*225, Ori_H*225);
disp([ 'PSNR:',num2str(psnr), ', SSIM:',num2str(ssim)]);

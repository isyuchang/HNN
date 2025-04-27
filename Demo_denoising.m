clear,clc

load('Pavia');


[M,N,B] = size(Ori_H);

sigma = 0.1; 
noise = sigma * randn(size(Ori_H)); 
Noi_H = Ori_H + noise;

[L2,S2] = HRPCA2D(Noi_H,'maxIter',200,'rho',1.05,'lambda',2,'debug',1,'GT',Ori_H,'tol', 1e-6);

[psnr, ssim, fsim, ergas, msam] = MSIQA(L2*225, Ori_H*225);
disp([ 'PSNR:',num2str(psnr), ', SSIM:',num2str(ssim)]);


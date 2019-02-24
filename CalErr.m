function [f]=CalErr(bhat0)
%Calculate error for the estimation bhat0

global Ividmeas HStack;

%bhat0=gpuArray(bhat0);

EstIntStack=abs(ifft2(bsxfun(@times,HStack,bhat0))).^2;%Fresnel propagation bhat0, inverse Fourier transform and take intensity
f=sum(sum(sum((Ividmeas-EstIntStack).^2)));




% tic
% for nz=1:Nz
%     
%     bhat=HStack(:,:,nz).*bhat0;
%     ahat= Ft(bhat);
%     EstInt=abs(ahat.^2);
%     diff(:,:,nz)=Ividmeas(:,:,nz)-EstInt;
%     
%     % diff(:,:,nz)=Ividmeas(:,:,nz)-(abs(Ft(HStack(:,:,nz).*bhat0)).^2);
% 
% end
% toc





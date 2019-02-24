function [y]=JHJx(x)
% Calculate y=(J^H)J*x


global  HStack ahatStack;

[Nx,Ny,Nz]=size(HStack);

N=Nx*Ny;


y=sum(conj(HStack).*fft2(ahatStack.*(real(conj(ahatStack).*(ifft2(bsxfun(@times,HStack,x)))))),3)/(N/2);



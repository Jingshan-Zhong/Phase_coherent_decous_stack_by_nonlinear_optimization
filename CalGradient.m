function [dfda]=CalGradient(bhat0)
%Gradient method for phase retrieval

global Ividmeas HStack ahatStack;

[Nx,Ny,Nz]=size(Ividmeas);

%bhat0=gpuArray(bhat0);

N=Nx*Ny;
%N=gpuArray(N);

ahatStack=ifft2(bsxfun(@times,HStack,bhat0));
ErrStack=Ividmeas-abs(ahatStack).^2;
dfdaStack=conj(HStack).*fft2(ahatStack.*ErrStack);
dfda=(-2/N)*sum(dfdaStack,3);



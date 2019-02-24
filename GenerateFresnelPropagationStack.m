function [HStack]=GenerateFresnelPropagationStack(Nx,Ny,z, nfocus, lambda, ps)


%=gpuArray();

cx=floor(Nx/2)+1;cy=floor(Ny/2)+1;
[us, vs]=ndgrid([1:Nx]-cx,[1:Ny]-cy);
us=us/Nx/ps; vs=vs/Ny/ps;
us=ifftshift(us); vs=ifftshift(vs);

Nz=length(z);

if isa(z,'gpuArray')
HStack=gpuArray.zeros(Nx,Ny,Nz);
else
HStack=zeros(Nx,Ny,Nz);
end

for nz=1:Nz
    
    HStack(:,:,nz)= exp(-1i*lambda*pi*(us.^2+vs.^2)*(z(nz)-z(nfocus)));

    
end


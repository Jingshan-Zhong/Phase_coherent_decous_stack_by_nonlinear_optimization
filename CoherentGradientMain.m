clear all; close all;
 
global Ividmeas HStack;
F = @(x) fft2(x);
Ft = @(x) ifft2(x);

 %% %%% Loading defocus data

load('SimulationCoherentDefocusStack');

IsGPUprogram=0; %choose whether run on gpu

%% Gradient method to solve phase
[Nx, Ny, Nz]=size(Ividmeas);

%%%%%%%%%%%%%load the parameters into GPU if IsGPUprogram=1 %%%%%%%%%%%%%%%%%%%%%%%%%%
if IsGPUprogram==1
Ividmeas = gpuArray(Ividmeas); Nx=gpuArray(Nx);Ny=gpuArray(Ny);z=gpuArray(z); nfocus=gpuArray(nfocus);lambda=gpuArray(lambda);ps=gpuArray(ps);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[HStack]=GenerateFresnelPropagationStack(Nx,Ny,z, nfocus, lambda, ps);

bhat0=F(Ividmeas(:,:,nfocus).^(1/2));%initialize the complex field with square root of intensity at focus
 
%parameters for line search
c=0.5;gamma=0.5;eps1=10^-4;

%%%%%%%%%%%%%%%%%%%%%%load the parameters into GPU if IsGPUprogram=1 %%%%%%%%%%%%%%%%%
if IsGPUprogram==1  
bhat0=gpuArray(bhat0);c=gpuArray(c);gamma=gpuArray(gamma);eps1=gpuArray(eps1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tic
MaxIter=20; maxiterCG=50;
bhat0=IterativeOptimization(@CalErr,@CalGradient,bhat0,MaxIter,c,gamma,eps1,maxiterCG);

%ltime=toc

ahat0=ifft2(bhat0);
ahat0=gather(ahat0);


%% show result

figure;
imagesc(Ividmeas(:,:,Nz));axis image;axis off;colormap gray
title(sprintf('Measured Intensity at most defocus distance'));colorbar

figure;
imagesc(abs(ahat0.^2));
axis image;axis off;colormap gray
title('Estimated intensity');colorbar

figure;
%imagesc(-angle(ahat0),[-0.5 1]);
imagesc(angle(ahat0));
axis image;axis off;colormap gray
title('Estimated phase');colorbar

% figure;
% imagesc((angle(ahat0)-mean(mean(angle(ahat0))))-(angle(TrueImg)-mean(mean(angle(TrueImg)))));colormap(gray);axis ij;axis equal;axis tight;
% %imagesc(angle(img)+angle(TrueImg(1:end-1,1:end-1)));colormap(gray);axis ij;axis equal;axis tight;
% axis image;axis off;colormap gray;colorbar;
% axis off;title('Phase Error');drawnow;
% 
% PhaseErr=angle(ahat0)-mean(mean(angle(ahat0)))-angle(TrueImg)+mean(mean(angle(TrueImg)));
% PHMSE=sum(sum(abs(PhaseErr)))



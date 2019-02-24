function bhat0=IterativeOptimization(f,df,bhat0,Ite,c,gamma,eps1,maxiterCG)
%the old name of this function is CoherentPhaseSteepGradient
%Gradient method for phase retrieval


if isa(bhat0,'gpuArray')
    Err=gpuArray.zeros(Ite,1);
    fcallstack=gpuArray.zeros(Ite,1);
    CGIte=gpuArray.zeros(Ite,1);
    Estbhat0Ite=gpuArray.zeros([size(bhat0) Ite]);
    
else
    Err=zeros(Ite,1);
    fcallstack=zeros(Ite,1);
    CGIte=zeros(Ite,1);
    Estbhat0Ite=zeros([size(bhat0) Ite]);
end

fk=feval(f,bhat0);
Err(1)=fk;
Estbhat0Ite(:,:,1)=ifft2(bhat0);

[Nx Ny]=size(bhat0); N=Nx*Ny;

%tic

GDinitialize=0;


for k=1:Ite
    
    [dfda]=feval(df,bhat0);
    dfdasum=sqrt(sum(sum(abs(dfda).^2)));
    %    if dfdasum <= StopCond
    %         disp('Meeting stopping condtions in step gradient method')
    %         break
    %     end
    

    if k<GDinitialize
        searchdirection=-dfda*N;%steep gradient  
         %Gauss Newton
        %
    else
        b=-0.5*dfda; CGinit=-dfda*N; %maxiterCG=50;
       
        [searchdirection, niter, flag] = solveCG(@JHJx, b, CGinit, norm(b)*10^-3, maxiterCG);
        CGIte(k)=niter;
        %}
    end
    
    %line search
    DDfnc=real(sum(sum(conj(searchdirection).*dfda)));
    display(DDfnc);
    
    if DDfnc>0
        searchdirection=b;
    end
    %}
    
    [bhat0,fk,fcall] = backtrack(bhat0,searchdirection,fk,f,DDfnc,c,gamma,eps1);
    fcallstack(k)=fcall;
    
    toc
    display(sprintf('Iteration= %d, Err=%f, fcallnumber=%d',k,fk,fcall));
    
    %save iteration results and number
    Err(k+1)=fk;
    %Estbhat0Ite(:,:,k+1)=ifft2(bhat0);
    
    
end


%save('-v7.3','CoherentNonOptError.mat','Err','Estbhat0Ite');

% figure;
% plot([1:length(Err)],Err);
%
% figure;
% plot([1:length(fcallstack)],fcallstack);
%
% figure;
% plot([1:length(CGIte)],CGIte);


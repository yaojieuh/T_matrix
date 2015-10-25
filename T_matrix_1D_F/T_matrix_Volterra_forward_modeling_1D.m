%% 1D acoustic wave modeling in T-matrix representation
%   the code have 2 choice: one for directly matrix inversion
%   the second is Volterra renormlization
%   Author: Jie Yao
%   Date:   10/23/2015


clc;clear;

% parpool('local',2)



% space grid point
Nx=401;Lx=1000;dx=Lx/(Nx-1);
xlable=[1:Nx]*dx-dx;
% wavelet initialization 
f0=25;
dt=0.004;Nt=301;
t0 = 1/f0;
lt=dt*(Nt-1);
[rt,time]=rickerwavelet(f0,dt,Nt,t0,1);
rw=ifft(rt);
Nw=(Nt-1)/2+1;
omega=2*pi/lt*[0:Nw-1 -Nw+1:-1];

%velocity initializaiton
c0=1500;
V=zeros(Nx,1);
for ix=1:Nx
    x=(ix-1)*dx;
    if (x>=300)&(x<600)
        V(ix)=0.5;
    end
end
Vmat=diag(V);

%Linear equation solver
solver='Matrix_inversion';  % 'Matrix_inversion' or 'Volterra'
 solver='Volterra';  % 'Matrix_inversion' or 'Volterra'

psir=zeros(Nx,Nt);

 tic
for iw=1:Nt
    k=omega(iw)/c0;
    for ix=1:Nx
          x=(ix-1)*dx;
          G1d(ix)=-1i*k/2*exp(1i*k*abs(x));
    end
    Gmat=toeplitz(G1d,G1d);
    switch lower(solver)
    case {'matrix_inversion'}
      
%       for ix=1:Nx        
%           for ix2=1:Nx
%               Gmat(ix,ix2)=-1i*k/2*exp(1i*k*abs(ix-ix2)*dx);
%           end
%       end   
      Imat= eye(Nx);
      A=Imat-Vmat*Gmat*dx;
      Tmat=inv(A)*Vmat;
      
    case {'volterra'}
        GmatV=zeros(Nx,Nx);
        for ix=1:Nx
          x=(ix-1)*dx;
          Gv1d(ix)=k*sin(k*x);
        end
        GmatV=tril(toeplitz(Gv1d));
%         for ix=1:Nx
%             for ix2=1:Nx
%                % Gmat(ix,ix2)=-1i*k/2*exp(1i*k*abs(ix-ix2)*dx);
%                 if(ix2<ix)              
%                     GmatV(ix,ix2)=k*sin(k*(ix-ix2)*dx);
%                 end      
%             end
%         end
       

        Imat= eye(Nx);
        A=Imat-Vmat*GmatV*dx;
%         invA=inv(A);
%         T1mat=invA*Vmat;
        T1mat=A\Vmat;
        
        for ix=1:Nx
            x=(ix-1)*dx;
            Gvect1(ix,1)=-1i*k/2*exp(-1i*k*x);
            Gvect2(1,ix)=exp(1i*k*x);
        end
        
%         Gvect3=Gvect1.';
%         T2vect=invA*(Vmat*Gvect3)*dx;
        T2vect=A\(Vmat*Gvect1*dx);
        B=(1-Gvect2*T2vect);
        Rvect=Gvect2*T1mat/B;
        Tmat=T1mat+T2vect*Rvect;
    otherwise
      disp(' solver undifined');
      return
    end
    
     for ix=1:Nx
        xr=(ix-1)*dx;    
        xs=0;
        psi0(ix)=-1i*k/2*exp(1i*k*abs(xs-xr));
        psis(ix)=-1i*k/2*exp(1i*k*(xr-xs));
        


     end
     psir(1:Nx,iw)=psi0.'+Gmat*Tmat*psis.'*dx;
     psir(1:Nx,iw)= psir(1:Nx,iw)*rw(iw);

end
toc
psirt=real(fft(psir,[],2));

figure
imagesc(time*1000,xlable,psirt)
ylabel('Distance (m)'); xlabel('Time (ms)');
title('Wave field');
axis image;colorbar
set(gca,'XAxisLocation','top');colormap(flipud(gray))
set(gca,'YTick',[0 500 1000]);
set(gca,'XTick',[0 500 1000]);



c=c0./sqrt(1-V);


%% inversion

% Nord=50;
% Va=zeros(Nx,Nx,Nord);
% Va(1:Nx,1:Nx,1)=Tmat;
% GV1=Gmat*Tmat;Vasum=Tmat;
% for it=2:Nord
%    Va(1:Nx,1:Nx,it)=-Va(1:Nx,1:Nx,it-1)*GV1*dx;
%    Vasum=Vasum+Va(1:Nx,1:Nx,it);
% end
% 
% Vapp1=diag(real(Vasum));
%  
% plot (Vapp1)
% 
% 
% 
% %scaling 
% epsion=0.5;
% Vas=zeros(Nx,Nx,Nord);
% Vas(1:Nx,1:Nx,1)=(1-epsion)*Tmat;
% GV1=Gmat*Vas(1:Nx,1:Nx,1);Vasums=Vas(1:Nx,1:Nx,1);
% for it=2:Nord
%    Vas(1:Nx,1:Nx,it)=-Vas(1:Nx,1:Nx,it-1)*GV1*dx+epsion*Vas(1:Nx,1:Nx,it-1);
%    Vasums=Vasums+Vas(1:Nx,1:Nx,it);
% end
% 
% 
% Vapp2=diag(real(Vasums));
% figure(2)
% plot (Vapp2)

% psi1(1)=exp(1i*k*0);psi1(2)=exp(1i*k*dx);
% for ix=3:Nx;
%     x=(ix-1)*dx;
%     psi1(ix)=exp(1i*k*x);
%     for ix2=1:ix-1;
%         x2=(ix2-1)*dx;
%         psi1(ix)=psi1(ix)+k*sin(k*(x-x2))*V(ix2)*psi1(ix2)*dx;
%     end
% end
% 
% I1=0;
% I2=0;
% for ix=1:Nx
%     x=(ix-1)*dx;
%     I1=I1+exp(1i*k*x)*V(ix)*psi1(ix)*dx;
%     I2=I2+exp(1i*k*x)*V(ix)*conj(psi1(ix))*dx;
% end
% Rk=-1i*k/2*I1/(1+1i*k/2*I2);
% psit=psi1+Rk*conj(psi1);
% 
% figure;
% plot (1:Nx,real(psir),1:Nx,real(psit),'r-')


% figure(2);
% plot (1:Nx,imag(psit),1:Nx,imag(psicheck),'ro')

% for ix=1:Nx
%     x=(ix-1)*dx;
%     psicheck(ix)=psis(ix);
%     psihalf(ix)=0;
%     for ix2=1:Nx
%         x1=(ix2-1)*dx;
% %         psicheck(ix)= psicheck(ix)-1i*k/2*exp(1i*k*abs(ix-ix2)*dx)*V(ix2)*psit(ix2)*dx;
%          psihalf(ix)=psihalf(ix)+Gmat(ix,ix2)*V(ix2)*psit(ix2)*dx;
%     end
%      psicheck(ix)= psicheck(ix)+psihalf(ix);
% end
% ooo=Gmat'*Vmat*psit'*dx;
% psicheck=psis+(Gmat*Vmat*psit.').'*dx;
% Vcheckl=Vmat*psit.';
% Vcheckr=Vmat*(psis+(Gmat*Vmat*psit.').'*dx).';
% 
% Tcheckl=Tmat*psis.';
% Tcheckr=(Vmat+Vmat*Gmat*Tmat*dx)*psis.';

% psicheck=psis+(Gmat*Vmat*psit.').'*dx;
% psicheck2=psis+(Gmat*Tmat*psis.').'*dx*dx;

% norm(Tmat-(Vmat+Vmat*Gmat*Tmat*dx))    


        
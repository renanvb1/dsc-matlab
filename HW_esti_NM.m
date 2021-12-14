% =========================================================================
%   Created by Renan Vieira on 01/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code  might be used as long as the author is properly credited.
%   This code aims to estimate layer permeabilities and
%   skin factors using the Nelder-Mead method.
% =========================================================================
clear; 
clc; close all     %limpando o console, a memoria e fechando os graficos
%%
% setting the unit conversion constants - DONT CHANGE!!!
global alphap alphat
alphap=19.03;
alphat=0.0003484;

% defining the initial time and the shut in time
t0=9.6e-7;
tp=960;
% defining the number of timesteps
dim=round(10*log10(tp/t0)+1);
% computing the time vector
t=fill_time_new(t0,tp,dim);
a=find(abs(t-9.6e-4)<1e-10);
% setting the wellbore radius
rw=0.0762;

% getting the selected case properties
cas='g';
[nlayers, phi, kx, kz, s, len, h]=get_data(cas);
% setting the wellbore flow-rate
qinj=1000;
% setting the oil viscosity
mio=3.1*ones(1,nlayers);
% setting the water viscosity
miw=0.52;%*ones(1,nlayers); 
% setting the total compressibility
ct=1.22e-4*ones(1,nlayers);

% by model hypothesis, ky will always be equal to kx
ky=kx;
% by hypothesis, rskin will always be equal to 0.56
rskin=0.56*ones(size(phi));
% computing kskin from S and rskin
kskin=sqrt(kx.*kz)./(s.*log(rskin/rw)+1);
% in all cases, no wellbore offcentering is considered
dz=h/2;
% computing the mean permeability
kj=(kx.*ky.*kz).^(1/3);
% getting the relative permeability data
for jj=1:nlayers
    % here, it will be considered a piston-like displacement
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(0);
end
% computing the hydraulic diffusivities
eta=kj.*kro(1)./phi./mio./ct;

% computing the endpoint mobilities
lohat=kro(1,:)./mio;
lwhat=krw(end,:)./miw;
% initializing the total mobility and fractional flow matrices
lt=zeros(length(kro(:,1)),nlayers);
dfw=lt;
% computing the total mobility and fractional flow derivative matrices
for jj=1:nlayers
    [lt(:,jj),dfw(:,jj)]=fill_data(length(sw(:,jj)),krw(:,jj),kro(:,jj),sw(:,jj),mio(jj),miw);
end

flagperiodo=1;
flap=1;
% computing the original pressure and flow-rate data
[ptrue,~,~,qjtrue]=HW_press_new(t,tp,flagperiodo,flap,kj,kx,ky,kz,h,...
    len,dz,kskin,rskin,s,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);

% discarding the data before t=9.6*e-3
pnoise=ptrue(a:dim);
qjnoise=qjtrue(a:dim,:);
taux=t(a:dim);

% reading the noise data files
erp=importdata('errpres.m');
erq=importdata('errqj.m');
% adjusting the size of the noise vectors
erp=erp(1:dim-a+1);
erq=erq(1:dim-a+1,1:nlayers);

% defining the sine amplitude
amp=1/100*max(pnoise);
% defining the sine period
per=12;
% defining the damping term;
dam=10;
% defining the flow-rate error standard deviation
stdq=0.01*max(qjnoise);

% adjusting the artificial noise
erp=amp*(sin(2*pi*taux/per)+erp/dam);
for jj=1:nlayers
    erq(:,jj)=erq(:,jj)*stdq(jj);
end
% adding the artificial noise to the true data
pnoise=pnoise+erp;
qjnoise=qjnoise+erq;
% computing the noisy pressure derivative
dpnoise=comp_deriv2(taux,pnoise,tp,1);

resmin=norm(pnoise-ptrue(a:dim))/norm(pnoise);
for jj=1:nlayers
    resmin=resmin+norm(qjnoise(:,jj)-qjtrue(a:dim,jj))/norm(qjnoise(:,jj));
end
fprintf('Residuo calculado usando os parametros corretos: %.3f \n', resmin);

[kxnm, kznm, snm, pnm, qjnm]=NelderMead(t,tp,flagperiodo,h,len,dz,kskin,...
    rskin,phi,qinj,dfw,ct,lt,lohat,lwhat,nlayers,rw,pnoise,qjnoise);
pnm=pnm(a:dim);
qjnm=qjnm(a:dim,:);
dpnm=comp_deriv2(taux,pnm,tp,1);

% plotting the noisy pressure and pressure derivative data
linw=1.5;
fonts=14;
x1=1e-3;
x2=10^(floor(log10(tp))+1);

% figure
% loglog(taux,pnoise,'vk',taux,dpnoise,'vk','LineWidth',linw)
% xlim([x1 x2])
% title('Reference Pressure Data (Including the Artificial Noise)', 'fontsize',fonts)

figure
loglog(taux,pnoise,'vk',taux,pnm,'-g','LineWidth',linw)
hold on
loglog(taux,dpnoise,'vk',taux,dpnm,'-g','LineWidth',linw)
xlim([x1 x2])
hold off
title('Ref Pressure vs NM Pressure', 'fontsize',fonts)

figure
semilogx(taux,qjnoise,'v', taux, qjnm,'-', 'LineWidth',linw)
xlim([x1 x2])
title('Ref Flow-Rate vs NM Flow-Rate', 'fontsize',fonts)


clear amp dam per stdq dim alphap alphat eta fonts linw lt t0 taux t x1 x2
clear jj ky kj dz dfw kskin rskin

%%
% function that applies the Nelder Mead method to estimate layer permeabilities and skin factors
function [kxest, kzest, Sjest, pest, qjest]=NelderMead(t,tp,flagperiodo,...
    h,len,dz,kskin,rskin,phi,qinj,dfw,ct,lt,lohat,lwhat,nlayers,rw,...
    pref,qref)
    
    dim=length(t);
    crop=dim+1-length(pref);
    % 1 para norma absoluta, 0 para norma relativa
    fobj=0;
    
    % defining the reflection, contraction and expansion coefficients
    alpha=1.0;
    beta=0.5;
    gamma=2.0;
    % defining the convergence acceptable tolerance
    maxdres=1e-15;%eps;%1e-10;%
    % defining the maximum acceptble iterations number
    imax=4e3;
    % setting the problem dimension (ie, the number of variables to be determined)
    n=3*nlayers;
    % initializing the simplex with zeros
    p=zeros(n+1,n);
    
    % defining the lower and upper boundaries for the estimated parameters
    lbound=zeros(n,1); 
    ubound=lbound;
    % considera o vetor como sendo kxs kzs Sjs
    for iaux=1:nlayers
        lbound(0*nlayers+iaux)=0.1;      % lower boundary equal to 0.1
        ubound(0*nlayers+iaux)=3e3;     % upper boundary equal to 3000
        lbound(1*nlayers+iaux)=0.1;    % lower boundary equal to 0.1
        ubound(1*nlayers+iaux)=3e3;   % upper boundary equal to 3000
        lbound(2*nlayers+iaux)=-1.0;  % lower boundary equal to -1
        ubound(2*nlayers+iaux)=12.0;   % upper boundary equal to 10
    end
       
    for iaux=1:n
        % filling the initial simplex with random numbers within the bound.
        p(:,iaux)=(ubound(iaux)-lbound(iaux))*rand(n+1,1)+lbound(iaux);
    end
    
    % computing the initial simplex function values
    y=zeros(n+1,1);
    for iaux=1:n+1
        [kx,kz,s,kj,eta]=init_point(p(iaux,:),lohat,phi,ct);
        [presaux,~,~,qjaux]=HW_press_new(t,tp,flagperiodo,kj,kx,kx,kz,...
            h,len,dz,kskin,rskin,s,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
        if fobj==1
            y(iaux)=norm(pref-presaux(crop:dim));
            for jj=1:nlayers
                y(iaux)=y(iaux)+norm(qref(:,jj)-qjaux(crop:dim,jj));
            end
        else
            y(iaux)=norm(pref-presaux(crop:dim))/norm(pref);
            for jj=1:nlayers
                y(iaux)=y(iaux)+norm(qref(:,jj)-qjaux(crop:dim,jj))/norm(qref(:,jj));
            end
        end
    end
    
    % initializing the iteration counter as 1
    icount=1;
    res=zeros(1,imax);
    % displaying the values of the worst initial simplex vertex
    ih=findmax(y);
    ph=p(ih,:);
    % initializing the residual with a high value
    res(1)=y(ih);
    dres=-res(1);
    fprintf(' Pior vertice do simplex inicial:\n kxj \t kvj \t Sj\n');
    for iaux=1:nlayers
        fprintf(' %.1f \t %.1f \t %.2f \n',ph(0*nlayers+iaux),ph(1*nlayers+iaux),ph(2*nlayers+iaux));
    end
    
    tol=0.1*res(1);
    
    % entering the iteration loop (based on fig 6 from Luersen)
%     while(icount<imax && res(max(icount,1))>tol && dres<-maxdres)
%     while(icount<imax && (res(max(icount-1,1))>tol || dres<-maxdres))
    while(icount<imax && dres<-maxdres)
        % finding ph (the point that implies in maximum y)
        ih=findmax(y);
        ph=p(ih,:);
        yh=y(ih);
        
        % finding pl (the point that implies in minimum y)
        il=findmin(y);
        pl=p(il,:);
        yl=y(il);
        
        % computing the centroid Pbar
        pbar=(sum(p)-ph)/(n);
        
        % performing the reflection operation
        pr=(1+alpha)*pbar-alpha*ph;
        % checking if the reflected point is within the boundaries
        for iaux=1:n
            if pr(iaux)<lbound(iaux)
                pr(iaux)=lbound(iaux);
            elseif pr(iaux)>ubound(iaux)
                pr(iaux)=ubound(iaux);
            end
        end
        % evaluating the function at pr
        [kx,kz,s,kj,eta]=init_point(pr,lohat,phi,ct);
        [presaux,~,~,qjaux]=HW_press_new(t,tp,flagperiodo,kj,kx,kx,kz,...
            h,len,dz,kskin,rskin,s,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
        if fobj==1
            yr=norm(pref-presaux(crop:dim));
            for jj=1:nlayers
                yr=yr+norm(qref(:,jj)-qjaux(crop:dim,jj));
            end
        else
            yr=norm(pref-presaux(crop:dim))/norm(pref);
            for jj=1:nlayers
                yr=yr+norm(qref(:,jj)-qjaux(crop:dim,jj))/norm(qref(:,jj));
            end
        end
        
        % entering the decision tree (fig. 6 from Luersen)
        if (yr<yl)
            % performing the expansion
            pe=gamma*pr+(1-gamma)*pbar;
            % checking if the reflected point is within the boundaries
            for iaux=1:n
                if pe(iaux)<lbound(iaux)
                    pe(iaux)=lbound(iaux);
                elseif pe(iaux)>ubound(iaux)
                    pe(iaux)=ubound(iaux);
                end
            end
            % evaluating the function at pe
            [kx,kz,s,kj,eta]=init_point(pe,lohat,phi,ct);
            [presaux,~,~,qjaux]=HW_press_new(t,tp,flagperiodo,kj,kx,kx,kz,...
                h,len,dz,kskin,rskin,s,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
            if fobj==1
                ye=norm(pref-presaux(crop:dim));
                for jj=1:nlayers
                    ye=ye+norm(qref(:,jj)-qjaux(crop:dim,jj));
                end
            else
                ye=norm(pref-presaux(crop:dim))/norm(pref);
                for jj=1:nlayers
                    ye=ye+norm(qref(:,jj)-qjaux(crop:dim,jj))/norm(qref(:,jj));
                end
            end
            % updating the simplex according to the expansion and reflection results
            if (ye<yl)
                p(ih,:)=pe;
                y(ih)=ye;
            else
                p(ih,:)=pr;
                y(ih)=yr;
            end
        else
%             % initializing cond as true before checking if yr is a new max
%             cond=true;
%             % checking if yr is a new maximum
%             for iaux=1:n+1
%                 % if yr is lower than any y(i), then setting cond as false
%                 if (iaux~=ih && y(iaux)>yr)
%                     cond=false;
%                     break
%                 end
%             end
            if y(ih)>yr
                cond=false;
            else
                cond=true;
            end
            % if yr is a new maximum, perform the contraction
            if cond
                % if the yr is lower than yh, replace ph by pr
                if yr<yh
                    ph=pr;
                    yh=yr;
                end
                % performing the contraction
                pc=beta*ph+(1-beta)*pbar;
                % checking if the contracted point is within the boundaries
                for iaux=1:n
                    if pc(iaux)<lbound(iaux)
                        pc(iaux)=lbound(iaux);
                    elseif pc(iaux)>ubound(iaux)
                        pc(iaux)=ubound(iaux);
                    end
                end
                % evaluating the function at pc
                [kx,kz,s,kj,eta]=init_point(pc,lohat,phi,ct);
                [presaux,~,~,qjaux]=HW_press_new(t,tp,flagperiodo,kj,kx,kx,kz,...
                    h,len,dz,kskin,rskin,s,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
                if fobj==1
                    yc=norm(pref-presaux(crop:dim));
                    for jj=1:nlayers
                        yc=yc+norm(qref(:,jj)-qjaux(crop:dim,jj));
                    end
                else
                    yc=norm(pref-presaux(crop:dim))/norm(pref);
                    for jj=1:nlayers
                        yc=yc+norm(qref(:,jj)-qjaux(crop:dim,jj))/norm(qref(:,jj));
                    end
                end
                if(yc>yh)
                    % for a failed contraction, replace all simplex points
                    for iaux=1:n+1
                        if iaux~=il
                            p(iaux,:)=(p(iaux,:)+pl)/2;
                        end
                    end
%                     ifail=ifail+1;
                else
                    % for a successful contraction, replace ph by pc
                    p(ih,:)=pc;
                    y(ih)=yc;
                end
            % if yr is not a new maximum, replace ph by pr
            else
                p(ih,:)=pr;
                y(ih)=yr;
            end
        end
        
        res(icount)=max(y);
        if icount-1>0
            dres=res(icount)-res(icount-1);
        end
        
        if (mod(icount,100)==0)
%             if (mod(log10(icount),1)==0)
            fprintf('Itera: %d Res: %.4f\n', icount, res(icount-1));
        end
        
        % updating the counter
        icount=icount+1;
        
    end
    % Printing a message to the user telling if the method converged or not
    if (icount==imax)
        fprintf('O método NAO convergiu :(\n');
        fprintf('Residuo final: %.5f \n', res(icount-1));
%     elseif res(icount-1)<tol
%         fprintf('O método convergiu :D \n');
%         fprintf(' No de iteracoes: %d \n',icount);
%         fprintf('Residuo final: %.6f \n', res(icount));
    else
        fprintf('Residuo persistente :| \n');
        fprintf(' No de iteracoes: %d \n',icount);
        fprintf('Residuo final: %.6f \n', res(icount-1));
        fprintf('Delta residuo: %.6f \n', dres);
    end
    
    % computing the pressure response using the estimated parameters
    [kxest,kzest,Sjest,~,~]=init_point(p(il,:),lohat,phi,ct);
    [pest,~,~,qjest]=HW_press_new(t,tp,flagperiodo,kj,kx,kx,kz,...
        h,len,dz,kskin,rskin,s,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
%     kzest=sqrt(kzest.*kxest);
    % plotting the residuals
    figure
    semilogy(res)
    
    fprintf(' Estimated parameters:\n kxj \t kzj \t Sj\n');
    for iaux=1:nlayers
        fprintf(' %.1f \t %.1f \t %.2f \n',kxest(iaux),kzest(iaux),Sjest(iaux));
    end
end

%%
% function that returns the index of the maximum entry of a vector
function [imax]=findmax(x)
    % initializing i as 1
    imax=1;
    if length(x)>1
        for i=2:length(x)
            if x(i)>x(imax)
                imax=i;
            end
        end
    end
end

%%
% function that returns the index of the maximum entry of a vector
function [imin]=findmin(x)
    % initializing i as 1
    imin=1;
    if length(x)>1
        for i=2:length(x)
            if x(i)<x(imin)
                imin=i;
            end
        end
    end
end

%%
% function gets as input a simplex point and returns some reservoir paramet
function [kx,kz,Sj,kj,etaj]=init_point(p,lohat,phi,ct)
    n=length(phi);
    kx=p(0*n+1:1*n);
    kz=p(1*n+1:2*n);
    Sj=p(2*n+1:3*n);
    kj=(kx.*kx.*kz).^(1/3.0);
    etaj=kj.*lohat./phi./ct;
end

%%
% function that gets the reservoir properties for each case
function [nlayers, phi, kx, kz, Sj, len, h]=get_data(flag)
    switch flag
        case 'a'
            nlayers=1;
            phi=0.23;
            kx=500;
            kz=500;
            Sj=0.0;
            len=700.0;
            h=20.0;
        case 'b'
            nlayers=1;
            phi=0.23;
            kx=500;
            kz=500;
            Sj=5.0;
            len=700.0;
            h=20.0;
        case 'c' % 2 cam 1 isotropica e 1 anisotropica
            nlayers=2;
            phi=[0.23 0.18];
            kx=[500 300];
            kz=[500 30];
            Sj=[5.0 4.0];
            len=[700 1000];
            h=[20 15];
        case 'd' %2 cam razao de anisotropia 1/1000
            nlayers=2;
            phi=[0.23 0.18];
            kx=[2000 1000];
            kz=[2 1];
            Sj=[5.0 9.0];
            len=[700 1000];
            h=[20 15];
        case 'e' %3 cam razoes de anisotropia 1/4 1/9 e 1/100
            nlayers=3;
            phi=[0.23 0.18 0.14];
            kx=[500 700 900];
            kz=[500 77.777 9.0];
            Sj=[5.0 7.0 6.0];
            len=[700 1000 1200];
            h=[20 15 10];
        case 'f' %2 cam razoes de anisotropia 1/400 e 1/1000 (redundante com o caso D)
            nlayers=2;
            phi=[0.23 0.18];
            kx=[800 1000];
            kz=[2 1];
            Sj=[5.0 7.0];
            len=[700 1000];
            h=[20 15];
        case 'g' %2 cam sem conseguir identificar o erf
            nlayers=2;
            phi=[0.23 0.18];
            kx=[800 1000];
            kz=[800 1000];
            Sj=[3.0 4.0];
            len=[200 200];
            h=[10 10];
        case 'h' %2 cam sem conseguir identificar o lrf
            nlayers=2;
            phi=[0.23 0.18];
            kx=[200 100];
            kz=[2 1];
            Sj=[5.0 7.0];
            len=[1000 1500];
            h=[20 15];
    end
end
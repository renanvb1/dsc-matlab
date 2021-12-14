% =========================================================================
%   Created by Renan Vieira on 02/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code  might be used as long as the author is properly credited.
%   This code plots the graphs required by the paper about estimating
%   individual layer properties from well test data in HW.
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
cas='e';
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
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(-1);
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

%%
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

% resmin=norm(pnoise-ptrue(a:dim))/norm(pnoise);
resmin=norm(pnoise-ptrue(a:dim));
for jj=1:nlayers
%     resmin=resmin+norm(qjnoise(:,jj)-qjtrue(a:dim,jj))/norm(qjnoise(:,jj));
    resmin=resmin+norm(qjnoise(:,jj)-qjtrue(a:dim,jj));
end
fprintf('Objective function computed using the true set of parameters: %.3f \n', resmin);

%% computing pressure and flow-rate profiles using the estimated properties
% getting the timesteps when ERF ends and when LRF starts
[terf, tlrf]=get_time(cas);

vas=1;
if vas==1
% getting the estimated properties
[kxrnpa, kzrnpa, srnpa, kxdtm, kzdtm, sdtm, kxnm, kznm, snm,...
    kzeq,kxeq,seq,steq]=get_esti(cas);

% computing the permeabilities in z direction
kzrnpa=kzrnpa.*kzrnpa./kxrnpa;
kzdtm=kzdtm.*kzdtm./kxdtm;
kznm=kznm.*kznm./kxnm;

% computing the average permeabilities
kjrnpa=(kxrnpa.*kxrnpa.*kzrnpa).^(1/3);
kjdtm=(kxdtm.*kxdtm.*kzdtm).^(1/3);
kjnm=(kxnm.*kxnm.*kznm).^(1/3);

% computing the mechanical skin
saux=sqrt(kxrnpa)+sqrt(kzrnpa);
saux=saux./sqrt(sqrt(kxrnpa.*kzrnpa));
saux=-log(saux/2);
srnpa=srnpa-saux;
saux=sqrt(kxdtm)+sqrt(kzdtm);
saux=saux./sqrt(sqrt(kxdtm.*kzdtm));
saux=-log(saux/2);
sdtm=sdtm-saux;

% computing the pressure and flow-rate data using rnpa
[prnpa,~,~,qjrnpa]=HW_press_new(t,tp,flagperiodo,flap,kjrnpa,kxrnpa,kxrnpa,kzrnpa,h,...
    len,dz,kskin,rskin,srnpa,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
% discarding the data before t=9.6*e-3
prnpa=prnpa(a:dim);
qjrnpa=qjrnpa(a:dim,:);

% computing the pressure and flow-rate data using dtm
[pdtm,~,~,qjdtm]=HW_press_new(t,tp,flagperiodo,flap,kjdtm,kxdtm,kxdtm,kzdtm,h,...
    len,dz,kskin,rskin,sdtm,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
% discarding the data before t=9.6*e-3
pdtm=pdtm(a:dim);
qjdtm=qjdtm(a:dim,:);

% computing the pressure and flow-rate data using nm
[pnm,~,~,qjnm]=HW_press_new(t,tp,flagperiodo,flap,kjnm,kxnm,kxnm,kznm,h,...
    len,dz,kskin,rskin,snm,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,1,1);
% discarding the data before t=9.6*e-3
pnm=pnm(a:dim);
qjnm=qjnm(a:dim,:);

% computing the derivatives
dprnpa=comp_deriv2(taux,prnpa,tp,1);
dpdtm=comp_deriv2(taux,pdtm,tp,1);
dpnm=comp_deriv2(taux,pnm,tp,1);
end
%% calculando as propriedades equivalentes do reservatorio
% permeabilidade equivalente no plano zx
kzx=sqrt(kz.*kx);
kzxeq=sum(len.*kzx)/sum(len);
% permeabilidade equivalente no plano xy
kxy=sqrt(kx.*ky);
kxyeq=sum(h.*kxy)/sum(h);

% mobilidade equivalente do oleo
lohatm=sum(lohat.*h)/sum(h);

% produto phi-ct equivalente
phicteq=sum(ct.*phi.*h)/sum(h);

% sani=-log((sqrt(kx)+sqrt(kz))./2./((kx.*kz).^0.25));
sani=sqrt(kx)+sqrt(kz);
sani=sani./sqrt(sqrt(kx.*kz));
sani=-log(sani/2)*1;
% Sj=Sj+sani;
% skin mecanico equivalente
seq=sum(kzx.*len.*(s+1*sani))/kzxeq/sum(len);

% expoente ej necessario para o calculo do skin total (lrf)
ej=kxy.*h./kzx./len;
% Sgj necessario para o calculo do skin total (lrf)
sgj=len.*((2*pi*rw./h).^ej);     sgj=log(4*rw./sgj);
% skin total em cada camada (que ocorre durante o lrf)
stj=sgj+kxy.*h./kzx./len.*s;
% skin total equivalente
steq=sum(kxy.*h.*stj)/kxyeq/sum(h);

% patamares teoricos de derivada
merf=alphap*qinj/2/lohatm/sum(len)/kzxeq;
mlrf=alphap*qinj/2/lohatm/sum(h)/kxyeq;

%% computing the early and late time logarithmic approximations
% getting the time range for ERF
t1=t(t<=terf);
% getting the time range for LRF
t2=t(t>=tlrf);

% definido o argumento do ln no curto tempo
argerf=4*kzxeq*lohatm*alphat/exp(.57722)/phicteq/rw/rw;
% definido o argumento do ln no longo tempo
arglrf=4*kxyeq*lohatm*alphat/exp(.57722)/phicteq/rw/rw;

% calculando as aproximacoes logaritmicas de curto e longo tempo
perf=merf*(log(argerf.*t1)+2*seq);
plrf=mlrf*(log(arglrf.*t2)+2*steq);


%% plotting the noisy pressure and pressure derivative data
linw=1.5;
fonts=14;
x1=1e-3;
x2=1000;
% y1=10^(floor(log10(min(pnoise)))-1);
y1=mean(dpnoise(1:15));
y1=10^(floor(log10(y1))-0);
y2=10^(ceil(log10(max(pnoise)))-0);
leg=strings(3,1);
% for jj=1:nlayers
%     leg(jj)=strcat('Layer ',' ',num2str(jj),' ');
% end
leg(1)='Layer 1';
leg(2)='Layer 2';
leg(3)='Layer 3';
leg(4)='Layer 4';
leg=leg(1:nlayers);

%% plotting the reference pressure data
figure
loglog(taux,pnoise,'vk','LineWidth',linw)
hold on
loglog(t1,perf,'*m','LineWidth',1.1*linw)
xline(terf,'--k','LineWidth',linw)
xline(tlrf,':k','LineWidth',linw)
loglog(t2,plrf,'*m','LineWidth',1.1*linw)
loglog(taux,dpnoise,'vk','LineWidth',linw)
hold off
xlim([x1 x2])
ylim([y1 y2])
xlabel('t (h)', 'fontsize',fonts)
ylabel('Pressure and press. deriv. (kgf/cm²)', 'fontsize',fonts)
legend({'Ref. Data', 'Log. Appr.', 'End of ERF', 'Start of LRF'}, 'Location','southeast', 'fontsize',fonts*.85)
% legend({'Ref. Data', 'Log. Appr.'}, 'Location','southeast', 'fontsize',fonts*.85)

%% plotting the reference flow-rate data
figure
semilogx(taux,qjnoise(:,1),'vk','LineWidth',linw)
hold on
semilogx(taux,qjnoise(:,2),'vb','LineWidth',linw)
if nlayers==3
    semilogx(taux,qjnoise(:,3),'vr','LineWidth',linw)
elseif nlayers ==4
    semilogx(taux,qjnoise(:,3),'vr','LineWidth',linw)
    semilogx(taux,qjnoise(:,4),'vg','LineWidth',linw)
end
hold off
legend(leg, 'Location','east', 'fontsize',fonts*.85)
xlabel('t (h)', 'fontsize',fonts)
ylabel('Flow-rate (m³/d)', 'fontsize',fonts)
xlim([x1 x2])

if vas==1
%% plotting the comparison between pressure data
figure
loglog(taux,pnoise,'vk','LineWidth',linw)
hold on
loglog(taux,prnpa,'xb','LineWidth',linw)
loglog(taux,pdtm,'dr','LineWidth',linw)
loglog(taux,pnm,'-g','LineWidth',linw)
loglog(taux,pnoise,'vk','LineWidth',linw)
loglog(taux,dprnpa,'xb','LineWidth',linw)
loglog(taux,dpdtm,'dr','LineWidth',linw)
loglog(taux,dpnm,'-g','LineWidth',linw)
loglog(taux,dpnoise,'vk','LineWidth',linw)
hold off
xlim([x1 x2])
ylim([y1 y2])
xlabel('t (h)', 'fontsize',fonts)
ylabel('Pressure and press. deriv. (kgf/cm²)', 'fontsize',fonts)
legend({'Ref.', 'RNPA', 'DTM', 'NM'}, 'Location','southeast', 'fontsize',fonts*.85)

%% plotting ref flow-rate data vs. nm flow-rate data
figure
semilogx(taux,qjnoise(:,1),'vk','LineWidth',linw)
hold on
semilogx(taux,qjnm(:,1),'-k','LineWidth',linw)
semilogx(taux,qjnoise(:,2),'vb','LineWidth',linw)
semilogx(taux,qjnm(:,2),'-b','LineWidth',linw)
if nlayers==3
    semilogx(taux,qjnoise(:,3),'vr','LineWidth',linw)
    semilogx(taux,qjnm(:,3),'-r','LineWidth',linw)
elseif nlayers==4
    semilogx(taux,qjnoise(:,3),'vr','LineWidth',linw)
    semilogx(taux,qjnm(:,3),'-r','LineWidth',linw)
    semilogx(taux,qjnoise(:,4),'vg','LineWidth',linw)
    semilogx(taux,qjnm(:,4),'-g','LineWidth',linw)
end
hold off
legend({'Ref.', 'NM'}, 'Location','east', 'fontsize',fonts*.85)
xlabel('t (h)', 'fontsize',fonts)
ylabel('Flow-rate (m³/d)', 'fontsize',fonts)
xlim([x1 x2])

%% plotting the comparison between flow-rate data
figure
semilogx(taux,qjnoise(:,1),'vk','LineWidth',linw)
hold on
semilogx(taux,qjrnpa(:,1),'xk','LineWidth',linw)
semilogx(taux,qjdtm(:,1),'dk','LineWidth',linw)
semilogx(taux,qjnm(:,1),'-k','LineWidth',linw)
semilogx(taux,qjnoise(:,2),'vb','LineWidth',linw)
semilogx(taux,qjrnpa(:,2),'xb','LineWidth',linw)
semilogx(taux,qjdtm(:,2),'db','LineWidth',linw)
semilogx(taux,qjnm(:,2),'-b','LineWidth',linw)
if nlayers==3
    semilogx(taux,qjnoise(:,3),'vr','LineWidth',linw)
    semilogx(taux,qjrnpa(:,3),'xr','LineWidth',linw)
    semilogx(taux,qjdtm(:,3),'dr','LineWidth',linw)
    semilogx(taux,qjnm(:,3),'-r','LineWidth',linw)
elseif nlayers==4
    semilogx(taux,qjnoise(:,3),'vr','LineWidth',linw)
    semilogx(taux,qjrnpa(:,3),'xr','LineWidth',linw)
    semilogx(taux,qjdtm(:,3),'dr','LineWidth',linw)
    semilogx(taux,qjnm(:,3),'-r','LineWidth',linw)
    semilogx(taux,qjnoise(:,4),'vg','LineWidth',linw)
    semilogx(taux,qjrnpa(:,4),'xg','LineWidth',linw)
    semilogx(taux,qjdtm(:,4),'dg','LineWidth',linw)
    semilogx(taux,qjnm(:,4),'-g','LineWidth',linw)
end
hold off
legend({'Ref.', 'RNPA', 'DTM', 'NM'}, 'Location','east', 'fontsize',fonts*.85)
xlabel('t (h)', 'fontsize',fonts)
ylabel('Flow-rate (m³/d)', 'fontsize',fonts)
xlim([x1 x2])
end
%%
clear alphap alphat amp ct dam dfw dz erp erq eta flagperiodo flap fonts jj
clear kj kjrnpa kjdtm kjnm kro krw kskin linw lt lwhat per x1 x2 y1 y2
clear a argerf arglrf dim ej leg lohatm miw perf plrf resmin rskin sani sgj
clear sw t1 t2 taux terf tlrf t0 vas
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
            kz=[125 77.777 9.0];
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
        case 'i' %3 cam alternativo ao caso e
            nlayers=3;
            phi=[0.23 0.18 0.14];
            kx=[500 800 900];
            kz=[62 50 9];
            Sj=[5 7 6];
            len=[700 1000 1000];
            h=[20 15 10];
        case 'j' %4 cam 
            nlayers=4;
            phi=[0.23 0.18 0.14 0.11];
            kx=[500 800 1000 600];
            kz=[56 50 1 6];
            Sj=[5 7 6 7];
            len=[600 1000 900 1000];
            h=[20 15 10 15];
    end
end

%%
% function that gets the moments when ERF ends and LRF starts for each case
function [terf, tlrf]=get_time(flag)
    switch flag
        case 'a'
            terf=.02411;
            tlrf=60.6;
        case 'b'
            terf=.02411;
            tlrf=60.6;
        case 'c' % 2 cam 1 isotropica e 1 anisotropica
            terf=.02411;
            tlrf=120.9;
        case 'd' %2 cam razao de anisotropia 1/1000
            terf=1.20857;
            tlrf=24.1;
        case 'e' %3 cam razoes de anisotropia 1/4 1/9 e 1/100
%             terf=.06057;
            terf=.02411;
            tlrf=76.3;
        case 'f' %2 cam razoes de anisotropia 1/400 e 1/1000 (redundante com o caso D)
            terf=1.20857;
            tlrf=76.3;
        case 'g' %2 cam sem conseguir identificar o erf
            terf=.00481;
            tlrf=2.41;
        case 'h' %2 cam sem conseguir identificar o lrf
            terf=1.20857;
            tlrf=481;
        case 'i' %3 cam alternativo ao caso e
            terf=.06057;
            tlrf=76.2;
        case 'j' %4 cam
            terf=.0763;
            tlrf=48.1;
    end
end

%%
% function that gets the estimated properties for each case (3 methods)
function [kxrnpa, kzrnpa, srnpa, kxdtm, kzdtm, sdtm, kxnm, kznm, snm, kzeq,kxeq,seq,steq]=get_esti(flag)
    switch flag
        case 'a'
            kxrnpa=549;
            kzrnpa=483;
            srnpa=-.6;
            kxdtm=531;
            kzdtm=648;
            sdtm=.2;
            kxnm=499;
            kznm=517;
            snm=.3;
            kzeq=481;
            kxeq=536;
            seq=-.62;
            steq=-7.41;
        case 'b'
            kxrnpa=550;
            kzrnpa=487;
            srnpa=4.3;
            kxdtm=532;
            kzdtm=8.8;
            sdtm=772;
            kxnm=499;
            kznm=517;
            snm=5.4;
            kzeq=480;
            kxeq=536;
            seq=4.17;
            steq=-7.26;
        case 'c' % 2 cam 1 isotropica e 1 anisotropica
            kxrnpa=[573 336];
            kzrnpa=[519 105];
            srnpa=[4.6 3.8];
            kxdtm=[519 321];
            kzdtm=[517 108];
            sdtm=[4.3 3.6];
            kxnm=[504 310];
            kznm=[559 99];
            snm=[6.8 4.7];
            kzeq=269;
            kxeq=455;
            seq=4.3;
            steq=-7.28;
        case 'd' %2 cam razao de anisotropia 1/1000
            kxrnpa=[2416 1079];
            kzrnpa=[70 32];
            srnpa=[3.7 7.5];
            kxdtm=[1996 1019];
            kzdtm=[66 30];
            sdtm=[3.1 6.5];
            kxnm=[2094 1086];
            kznm=[61 29];
            snm=[3.7 7.1];
            kzeq=46.3;
            kxeq=1729;
            seq=4.74;
            steq=.01;
        case 'e' %3 cam razoes de anisotropia 1/4 1/9 e 1/100
            vas=1; % 2 = std 1 = new, 0 = old
            if vas==0
            kxrnpa=[542 718 1001];
            kzrnpa=[255 236 80];
            srnpa=[4.5 6.2 4.0];
            kxdtm=[527 741 953];
            kzdtm=[265 249 89];
            sdtm=[4.5 6.5 4.5];
            kxnm=[507 723 960];
            kznm=[519 307 100];
            snm=[5.9 10.9 7.3];
            elseif vas==1
            kxrnpa=[542 716 1000];
            kzrnpa=[255 236 75];
            srnpa=[4.6 6.5 3.5];
            kxdtm=[526 740 952];
            kzdtm=[262 206 97];
            sdtm=[4.6 4.9 5.2];
            kxnm=[515 695 946];
            kznm=[255 230 86];
            snm=[5.0 6.7 4.8];
            else
            kxrnpa=[542 716 1000];
            kzrnpa=[255 236 75];
            srnpa=[4.6 6.5 3.5];
            kxdtm=[526 740 952];
            kzdtm=[262 206 97];
            sdtm=[4.6 4.9 5.2];
            kxnm=[501 704 922];
            kznm=[235 210 80];
            snm=[4.6 6.1 5.2];
            end
            kzeq=117.3;
            kxeq=687;
            seq=5.16;
            steq=-7.19;
        case 'f' %2 cam razoes de anisotropia 1/400 e 1/1000 (redundante com o caso D)
            kxrnpa=[927 1003];
            kzrnpa=[44 34];
            srnpa=[4.0 5.7];
            kxdtm=[889 1119];
            kzdtm=[43 34];
            sdtm=[3.5 5.4];
            kxnm=[811 1045];
            kznm=[39 31];
            snm=[4.1 5.6];
            kzeq=37;
            kxeq=948;
            seq=4.66;
            steq=-2.46;
        case 'g' %2 cam sem conseguir identificar o erf
            kxrnpa=[820 983];
            kzrnpa=[685 766];
            srnpa=[1.8 1.9];
            kxdtm=[801 1015];
            kzdtm=[646 689];
            sdtm=[1.2 1.1];
            kxnm=[794 986];
            kznm=[845 1097];
            snm=[3.4 4.7];
            kzeq=719;
            kxeq=901;
            seq=1.81;
            steq=-5.99;
        case 'h' %2 cam sem conseguir identificar o lrf
            kxrnpa=[374 121];
            kzrnpa=[21 10];
            srnpa=[4.1 6.0];
            kxdtm=[216 123];
            kzdtm=[21 10];
            sdtm=[4.0 5.6];
            kxnm=[201 102];
            kznm=[20 10];
            snm=[4.4 6.2];
            kzeq=14;
            kxeq=196;
            seq=4.64;
            steq=-6.1;
        case 'i' %3 cam alternativo ao caso e
            kxrnpa=[541 811 985];
            kzrnpa=[180 200 83];
            srnpa=[4.5 6.2 4.2];
            kxdtm=[526 841 939];
            kzdtm=[184 226 90];
            sdtm=[4.4 7.1 4.6];
            kxnm=[487 783 890];
            kznm=[183 203 87];
            snm=[5.0 6.7 5.1];
            kzeq=150;
            kxeq=715;
            seq=5.12;
            steq=-7.00;
        case 'j' %4 cam 
            kxrnpa=[524 798 1125 616];
            kzrnpa=[171 199 29 59];
            srnpa=[4.5 6.2 3.8 5.9];
            kxdtm=[514 830 1045 624];
            kzdtm=[174 225 32 69];
            sdtm=[4.4 7.0 4.3 6.9];
            kxnm=[482 778 964 572];
            kznm=[158 179 28 55];
            snm=[4.0 5.2 3.8 5.2];
            kzeq=108;
            kxeq=653;
            seq=5.29;
            steq=-6.58;
    end
end
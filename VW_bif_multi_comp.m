% =========================================================================
%   Created by Renan Vieira on 04/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   O codigo busca comparar o modelo desenvolvido usando funcoes
%   de Green para Tinj em VW multicamadas com a solucao existente, 
%   desenvolvida por Barreto et al. e Peres et al.
% =========================================================================

%% limpando o console e a memoria, e fechando as janelas graficas =========
clear; clc; close all     
tic
%% sessao de definicao do vetor de tempo e das variaveis globais ==========
% constante de conversao de unidades de pressao
alphap=19.03;
% constante de conversao de unidades de tempo
alphat=0.0003484;

% definido o flag que indica a existencia ou nao da estatica
fper=1;

% definindo o flag da CCE (0 = IAR; 1 = CPB; 2 = NFB; 3 = Canal; 4 = Falha)
fcce=0;

% definindo o tempo inicial
t0=4.8e-8;
% definindo o tempo de fechamento do poco (definido pelo usuario)
tp=9.6*1e+1;
% definindo o numero de passos em cada periodo 
dim=round(log10(tp/t0));
dim=5*dim+1;

%% calculando o vetor de tempo ============================================
t=fill_time(t0,tp,dim);
% se houver periodo de falloff, criar os passos de tempo para o falloff
if fper==2
    % definido o tempo total do teste 
    tend=2*tp;
    % chamando a funcao para criar um vetor auxiliar com os passos de tempo do falloff
    taux=fill_time(t0,tend-tp,dim);
    taux=taux+tp;
    % acrescentando os passos de tempo do falloff no vetor de tempo
    t=[t; taux];
    % limpando a variavel com o vetor de tempo auxiliar
    clear taux
end

%% lendo os dados de entrada no arquivo txt ===============================
filename='props2.txt';
fid=importdata(filename);
fid=fid.data;

% determinando o numero de camadas a partir do arquivo txt 
nlay=length(fid(:,1));
% lendo o raio do poco
rw=fid(1,1);
% lendo a vazao
q=fid(1,2);

% lendo as permeabilidades na regiao 1
ks=fid(:,3)';
% lendo os raios da regiao 1
rs=fid(:,4)';
% lendo as permeabilidades na regiao 2
k=fid(:,5)';
% lendo os raios da regiao 2
re=fid(:,6)';
% lendo as permeabilidades na direcao z
% kz=fid(:,7)';
% lendo as porosidades das camadas
phi=fid(:,8)';
% lendo as espessuras das camadas
h=fid(:,9)';
% lendo o comprimento do poco em cada camada
len=fid(:,10)';
% lendo as viscosidades do oleo
muo=fid(:,11)';

%% definindo alguns parametros de entrada =================================
% definindo o valor da viscosidade da agua (em cP)
muw=0.5;
% definindo o fator volume de formacao da agua (em m^3/ STD m^3)
bw=1;
% buscando os flags de permeabilidade relativa (definidos pelo usuario na interface)
flag=fid(:,12)';

%% calculando parametros que sao definidos a partir dos dados de entrada ==

for jj=1:nlay
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(flag(jj));
end
dfp=zeros(length(kro(:,1)),nlay);
ltp=dfp;
for jj=1:nlay
    [ltp(:,jj),dfp(:,jj)]=fill_data(length(sw(:,jj)),krw(:,jj),kro(:,jj),sw(:,jj),muo(jj),muw);
end
dfg=1./(sw(end,:)-sw(1,:));

% calculando as mobilidades extremais
logr=kro(1).*k./muo;
lwgr=krw(end).*ks./muw;

lop=logr./k;
lwp=lwgr./ks;

% calculando a compressibilidade total (ct = cr + cw*swi + co*(1-sor))
ct=2.9e-5*ones(size(h)); % cr+sw(1).*cw+sw(end).*co; %
% calculando a difusividade hidraulica (eta = k2*lohat/(phi*ct))
etao=alphat.*logr./phi./ct;
etaw=etao./logr.*lwgr;
% ajustando as unidades da vazao
q=alphap*q*bw;

%% sessao onde sera chamada a funcao que calcula a pressao e as vazoes ====
% calculando os dados de pressao e vazao
[pp, qp, rfp]=VW_peres(t,nlay,q,rw,dfp,sw,ltp,ct,phi,h,ks,k,rs,etaw,etao,lwp,lop);
[pg, qg, rfg]=VW_gre6(t,nlay,q,rw,dfg,phi,h,ks,k,rs,re,etaw,etao,lwgr,logr,fcce);

% calculando a derivada da pressao
[dpg]=compute_derivative(t,pg,tp,fper);
% calculando a derivada da pressao
[dpa]=compute_derivative(t,pp,tp,fper);

loeq=sum(logr.*h);
lweq=sum(lwgr.*h);

% calculando o patamar teorico da derivada 
mo=q/2/loeq;
mw=mo*loeq/lweq*ks/k;
% mw=qinj*miw/2/k/krw(end)/h;

q=q/alphap/bw;
clc;
%% plotando os resultados =================================================

% definindo a posicao da legenda
switch fcce
    case 0
        loc='east';
    case 1
        loc='east';
    case 2
        loc='northwest';
    case 3
        loc='northwest';
    case 4
        loc='northwest';
    case 5
        loc='east';
    otherwise
        loc='northwest';
end
% definindo as entradas da legenda
leg={'\Delta P (this work) ',  'der. (this work) ',...
    '\Delta P (Barreto et al.) ', 'der. (Barreto et al.) '};

% definindo os limites dos eixos nos graficos
bo=round(log10(t(1))); bo=max([10^(bo) 1e-4]); 
bw=floor(log10(tp))+1; bw=10^bw;

co=min(mo, mw); co=1*10^floor(log10(co));
cw=ceil(log10(max(pp))); cw=10^cw;
% co=cw/100;

linw=1.5;
fonts=14;
% plotando os dados durante o periodo de injecao
figure
% loglog(t,pwf,'-k','LineWidth',linw);
loglog(t,pg,'ok', t,dpg,'ob','LineWidth',linw);
hold on
loglog(t,pp,'-k', t,dpa,'-b','LineWidth',linw);
yline(mo,'--k')
yline(mw,':k','LineWidth',2*linw)
legend(leg,'Location',loc, 'fontsize',fonts*.75)
% legend({'\Delta P (canal) ', 'der. (canal) ','\Delta P (no-flow) ', ...
%     'der. (no-flow) '},'Location','northwest', 'fontsize',fonts*.75)
axis([bo,bw, co,cw])
title('Pressure and Pressure Derivative Profile')
xlabel('t (h)', 'fontsize',fonts)
ylabel('\Delta P (kgf/cm??)', 'fontsize',fonts)
grid on

figure; plot(qp)
hold on
plot(qg,'o')
plot(sum(qg'))
yline(q*alphap)

%% limpando algumas variaveis que nao sao mais necessarias ================
clear alphap alphat bo bw co cr cw dim dfw fid filename flag fper jj kro
clear krw p0 linw fonts loc etao etaw leg
% clear 
clc; toc

%% funcao que cria o vetor de tempo =======================================
function [t]=fill_time(t0,tend,npassos)
    % os inputs da funcao sao o no de passos, os tempos inicial e final    
    % inicializando o vetor de tempo com zeros
    t=zeros(npassos,1);
    % a funcao preenche o vetor de tempo conforme uma progressao geometrica
    tim=(tend/t0)^(1/(npassos-1));
    
    % definindo o tempo inicial na primeira posicao do vetor
    t(1)=t0;
    % preenchendo as outras posicoes do vetor conforme a PG
    for ii=2:npassos
        t(ii)=t(ii-1)*tim;
    end
end

%% nova funcao que calcula a pressao a cada passo de tempo ================
function [pwf,qj,rf]=VW_gre6(t,nlay,q,rw,dfw,phi,h,ks,k,rs,re,etaw,etao,lw,lo,fcce)
    % nao e necessario calcular as pressoes para tempos extremamente curtos
    % (apesar de a convolucao precisar que o tempo inicial seja pequeno)
    iaux=find(t>1e-5); iaux=iaux(1);
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay);
    pw2=pwf;
    
    % definindo o tamanho do infinito (usado para CCE canal)
    if fcce==3
        infi=50;
    end
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % por enquanto, vamos chutar as vazoes como sendo constantes
    for jj=1:nlay
        qj(:,jj)=q*k(jj)*h(jj)/sum(k.*h)*ones(dim,1);
    end
    
    fac=1.001;
    % inicializando os rFs em cada camada com zeros
    rf=zeros(1,nlay);
    % para cada tempo, calcular as respectivas fontes instantaneas
    for ii=iaux:dim
        % aproximando a vazao no passo de tempo atual pela vazao anterior
        qj(ii,:)=qj(ii-1,:);
        % calculando o raio da frente de avanco em cada camada
        rf=calc_rf(ii, rw, t, qj, h, phi, dfw);
        % inicializando algumas variaveis auxiliares
        r1=rf;          r2=r1;
        eta2=etaw;      l2=lw;
        % verificando se o rF ultrapassou o rskin
        for j=1:nlay
            if rf(j) < rs(j)
                r1(j)=rf(j); r2(j)=rs(j); 
                eta2(j)=etao(j)./k(j).*ks(j);   l2(j)=lo(j)./k(j).*ks(j);
            else
                r1(j)=rs(j); r2(j)=rf(j);
                eta2(j)=etaw(j).*k(j)./ks(j);   l2(j)=lw(j)./ks(j).*k(j);
            end
        end
        % entrando no loop da convolucao
        for ia=2:ii
            % inicializando os impulsos auxiliares com zeros
            g1=zeros(1,nlay);
            g2=zeros(1,nlay);
            g3=zeros(1,nlay);
            ga1=g1;
            ga2=g1;
            ga3=g1;
            % entrando no loop das camadas
            for jj=1:nlay
                % entrando no loop de stehfest
                for kk=1:N
                    % calculando a "variavel de Laplace"
                    u1=log(2)*kk/(t(ia-1));
                    u2=log(2)*kk/(t(ia+0)+t(ia-1))*2;
                    u3=log(2)*kk/(t(ia+0));
                    % calculando o valor de x
                    x1=Vj(N,kk)*sisgre(u1,1,h(jj),rw,r1(jj),r2(jj),re(jj),lw(jj),l2(jj),lo(jj),etaw(jj),eta2(jj),etao(jj),0,fcce);
                    x2=Vj(N,kk)*sisgre(u2,1,h(jj),rw,r1(jj),r2(jj),re(jj),lw(jj),l2(jj),lo(jj),etaw(jj),eta2(jj),etao(jj),0,fcce);
                    x3=Vj(N,kk)*sisgre(u3,1,h(jj),rw,r1(jj),r2(jj),re(jj),lw(jj),l2(jj),lo(jj),etaw(jj),eta2(jj),etao(jj),0,fcce);
                    % incrementando a soma nos impulsos auxiliares
                    g1(jj)=g1(jj)+x1(1)*besseli(0,sqrt(u1/etaw(jj))*rw);
                    g1(jj)=g1(jj)+x1(2)*besselk(0,sqrt(u1/etaw(jj))*rw);
                    
                    g2(jj)=g2(jj)+x2(1)*besseli(0,sqrt(u2/etaw(jj))*rw);
                    g2(jj)=g2(jj)+x2(2)*besselk(0,sqrt(u2/etaw(jj))*rw);
                    
                    g3(jj)=g3(jj)+x3(1)*besseli(0,sqrt(u3/etaw(jj))*rw);
                    g3(jj)=g3(jj)+x3(2)*besselk(0,sqrt(u3/etaw(jj))*rw);
                    
                    ga1(jj)=ga1(jj)+x1(1)*besseli(0,sqrt(u1/etaw(jj))*rw*fac);
                    ga1(jj)=ga1(jj)+x1(2)*besselk(0,sqrt(u1/etaw(jj))*rw*fac);
                    
                    ga2(jj)=ga2(jj)+x2(1)*besseli(0,sqrt(u2/etaw(jj))*rw*fac);
                    ga2(jj)=ga2(jj)+x2(2)*besselk(0,sqrt(u2/etaw(jj))*rw*fac);
                    
                    ga3(jj)=ga3(jj)+x3(1)*besseli(0,sqrt(u3/etaw(jj))*rw*fac);
                    ga3(jj)=ga3(jj)+x3(2)*besselk(0,sqrt(u3/etaw(jj))*rw*fac);
                    
                    switch fcce
                        % somando os pocos imagem para canal
                        case 3
                            for n=1:infi
                                g1(jj)=g1(jj)+2*x1(5)*besselk(0,sqrt(u1/etao(jj))*n*re(jj)*2);
                                g2(jj)=g2(jj)+2*x2(5)*besselk(0,sqrt(u2/etao(jj))*n*re(jj)*2);
                                g3(jj)=g3(jj)+2*x3(5)*besselk(0,sqrt(u3/etao(jj))*n*re(jj)*2);
                                
                                ga1(jj)=ga1(jj)+2*x1(5)*besselk(0,sqrt(u1/etao(jj))*n*re(jj)*2);
                                ga2(jj)=ga2(jj)+2*x2(5)*besselk(0,sqrt(u2/etao(jj))*n*re(jj)*2);
                                ga3(jj)=ga3(jj)+2*x3(5)*besselk(0,sqrt(u3/etao(jj))*n*re(jj)*2);
                            end
                            % somando o poco imagem para falha selante
                        case 4
                            g1(jj)=g1(jj)+x1(5)*besselk(0,sqrt(u1/etao(jj))*re(jj)*2);
                            g2(jj)=g2(jj)+x2(5)*besselk(0,sqrt(u2/etao(jj))*re(jj)*2);
                            g3(jj)=g3(jj)+x3(5)*besselk(0,sqrt(u3/etao(jj))*re(jj)*2);
                            
                            ga1(jj)=ga1(jj)+x1(5)*besselk(0,sqrt(u1/etao(jj))*re(jj)*2);
                            ga2(jj)=ga2(jj)+x2(5)*besselk(0,sqrt(u2/etao(jj))*re(jj)*2);
                            ga3(jj)=ga3(jj)+x3(5)*besselk(0,sqrt(u3/etao(jj))*re(jj)*2);
                            % somando o poco imagem para aquifero lateral
                        case 5
                            g1(jj)=g1(jj)-x1(5)*besselk(0,sqrt(u1/etao(jj))*re(jj)*2);
                            g2(jj)=g2(jj)-x2(5)*besselk(0,sqrt(u2/etao(jj))*re(jj)*2);
                            g3(jj)=g3(jj)-x3(5)*besselk(0,sqrt(u3/etao(jj))*re(jj)*2);
                            
                            ga1(jj)=ga1(jj)-x1(5)*besselk(0,sqrt(u1/etao(jj))*re(jj)*2);
                            ga2(jj)=ga2(jj)-x2(5)*besselk(0,sqrt(u2/etao(jj))*re(jj)*2);
                            ga3(jj)=ga3(jj)-x3(5)*besselk(0,sqrt(u3/etao(jj))*re(jj)*2);
                    end
                end
            end
            % calculando os integrandos com base nos impulsos auxiliares
            g1=1/sum(1./g1);
            g2=1/sum(1./g2);
            g3=1/sum(1./g3);
            
            ga1=1/sum(1./ga1);
            ga2=1/sum(1./ga2);
            ga3=1/sum(1./ga3);
            % forcando os termos NaN a ser iguais a zero
%             if isnan(g1)
%                 g1=0;
%             end
%             if isnan(g2)
%                 g2=0;
%             end
%             if isnan(g3)
%                 g3=0;
%             end
%             
%             if isnan(ga1)
%                 ga1=0;
%             end
%             if isnan(ga2)
%                 ga2=0;
%             end
%             if isnan(ga3)
%                 ga3=0;
%             end
            % invertendo os termos para o campo real
            g1=g1*log(2)/(t(ia-1));
            g2=g2*log(2)/(t(ia+0)+t(ia-1))*2;
            g3=g3*log(2)/(t(ia+0));
            
            ga1=ga1*log(2)/(t(ia-1));
            ga2=ga2*log(2)/(t(ia+0)+t(ia-1))*2;
            ga3=ga3*log(2)/(t(ia+0));
            
            % incrementando a integral
            pwf(ii)=pwf(ii)+(t(ia)-t(ia-1))*(g1+4*g2+g3)/6;
            pw2(ii)=pw2(ii)+(t(ia)-t(ia-1))*(ga1+4*ga2+ga3)/6;
        end
        % atualizando as vazoes usando uma derivada numerica
%         qj(ii,:)=(pwf(ii)-pw2(ii))/(fac-1).*lw.*h;
    end
    pwf=pwf*q;
    
end

%% sistema linear do impulso instantaneo radial no espaco de Laplace ======
function [x]=sisgre(u,Q,h,rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,t0,fcce)
    
    % definindo os argumentos das fun????es de bessel
    arg1=sqrt(rw*rw*u/eta1);
    arg2=sqrt(r1*r1*u/eta1);
    arg3=sqrt(r1*r1*u/eta2);
    arg4=sqrt(r2*r2*u/eta2);
    arg5=sqrt(r2*r2*u/eta3);
    
    c=l2/l1*sqrt(eta1/eta2);
    d=l3/l2*sqrt(eta2/eta3);
    
    % selecionando o sistema linear adequado com a CCE
    switch fcce    
        case 1    % sistema linear para CPB
            arg6=sqrt(re*re*u/eta3);
            % inicializando o vetor com zeros
            b=zeros(6,1);
            b(6)=1*0/u;
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0   0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)  0   0
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)  0   0
                0                0                  besseli(0,arg4)    besselk(0,arg4)   -besseli(0,arg5)   -besselk(0,arg5)
                0                0                  besseli(1,arg4)   -besselk(1,arg4) -d*besseli(1,arg5)  d*besselk(1,arg5)
                0                0                  0                  0                  besseli(0,arg6)    besselk(0,arg6)];
        
        case 2    % sistema linear para NFB
            arg6=sqrt(re*re*u/eta3);
            % inicializando o vetor com zeros
            b=zeros(6,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0   0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)  0   0
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)  0   0
                0                0                  besseli(0,arg4)    besselk(0,arg4)   -besseli(0,arg5)   -besselk(0,arg5)
                0                0                  besseli(1,arg4)   -besselk(1,arg4) -d*besseli(1,arg5)  d*besselk(1,arg5)
                0                0                  0                  0                  besseli(1,arg6)   -besselk(1,arg6)];
        
        otherwise    % sistema linear para IAR (tambem usado na CCE canal)
            % inicializando o vetor com zeros
            b=zeros(5,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)  0
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)  0
                0                0                  besseli(0,arg4)    besselk(0,arg4)  -besselk(0,arg5)
                0                0                  besseli(1,arg4)   -besselk(1,arg4) d*besselk(1,arg5)];
    end
    % alterando a 1a entrada do vetor
    b(1)=-Q/h/l1/arg1*exp(-u*t0);
    
    a=max(max(m));
    m=m/a;
    b=b/a;
    % calculando os coeficientes
    x=m\b;
%     x=inv(m)*b;
   
   % calculando a pressao no poco
%    y=x(1)*besseli(0,arg1)+x(2)*besselk(0,arg1);
   
end

%% funcao que calcula a pressao pelo modelo rad comp ======================
function [pwf, qj, rf]=VW_peres(t,nlay,q,rw,dfw,sw,lt,~,phi,h,ks,k,rs,~,etao,~,lo)
%     iaux=find(t>1e-5); iaux=iaux(1);
    iaux=2;
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de coeficientes Aj
    pwf = zeros(dim,1);
    qj = zeros(dim,1);
    aj = ones(1,nlay);
    rf = zeros(size(sw));
    
    % calculando os fatores de skin em cada camada
    s=(k./ks-1).*log(rs./rw);
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % por enquanto, vamos chutar as vazoes como sendo constantes
    for j=1:nlay
        qj(:,j)=q*k(j)*h(j)/sum(k.*h)*ones(dim,1);
    end
    
    % calculando a altura total do reservatorio
    ht=sum(h);
    % calculando a permeabilidade equivalente
    keq=k.*h/ht;
    % produto porosidade-compressibilidade equivalente
%     phict=phi.*ct.*h./ht;
    
    % entrando no loop de tempo
    for ii=iaux:dim
        % atualizando provisoriamente as vazoes com o valor do t anterior
        qj(ii,:)=qj(ii-1,:);
        
        % calculando o skin total do reservatorio
        seq=s.*qj(ii,:)/q;
        % calculando o termo deltaPskin
        dps=q*seq/sum(h.*lo)/keq;
        
        % calculando o termo deltapo
        dpo=0;
        % entrando no loop de stehfest
        for kk=1:N
            % calculando a "variavel de Laplace"
            u=log(2)*kk/t(ii);
            % calculando o valor de x
            x=Vj(N,kk)*lefk(u,rw,q,h,k,etao,lo);
            % incrementando a soma
            dpo=dpo+x;
        end
        % invertendo a soma para o campo real
        dpo=dpo*log(2)/t(ii);
        
        % calculando o termo deltaPlambda
        dpl=0;
        % entrando no loop das camadas
        for j=1:nlay
            % calculando o rF em cada camada
            rf(:,j)=calc_rf(ii, rw, t, qj(:,j), h(j), phi(j), dfw(:,j));
            % calculando o coeficiente aj para a camada
            aj(j)=calc_aj(rw,k(j),h(j),ks(j),rs(j),rf(:,j),lo(j),lt(:,j));
            % incrementando a soma
            dpl=dpl+1/aj(j);
        end
        % invertendo a soma
        dpl=q/dpl;
        
        % calculando a pressao no tempo atual
        pwf(ii)=dpo+dpl+dps;
        
        % atualizando as vazoes
        qj(ii,:)=(pwf(ii)-dpo-dps)./aj;
    end
    
end

%% funcao que calcula o coef Aj em um determinado tempo ===================
function [aj]=calc_aj(rw,k,h,ks,rs,rf,lo,lt)
    % determinando o comprimento do vetor sw
    comps=length(lt);
    % inicializando o aj
    aj=0.0;
    % calculando a integral necessaria para calcular aj
    for u=comps:-1:2
        % numerically integrating expression (lambdahat/ltotal - 1) d ln(r)/r using the trapeze rule
        aj=aj+log(rf(u-1)/rf(u))*(lo/lt(u)-1+lo/lt(u-1)-1)/2;
    end
    % if there is formation damage, applying the proper adjustment
    if rs>rw
        % if the waterfront is within the damaged zone, use equation 44 @pdf 29
        if rf(1) < rs
            % equation 44 @ pdf 29 from file "03 ... n-camadas.pdf"
            aj=aj/ks/h/lo;
        else
            % initializing the second integral in equation 45 @pdf 29
            intsj=0.0;
            % initializing u as the water saturation vector length
            u=comps;
            % incrementing the auxiliary integral for all radii smaller than rskinj
            while(rf(u)<=rs)
                % checking if the next radial step is out of the damaged zone
                if (rf(u-1)>rs)
                    % if the next radius is out of the damaged zone, interpolate the total mobility
                    ltskin=lt(u)+(lt(u-1) - lt(u))*(rs - rf(u)) / (rf(u-1) - rf(u));
                    % incrementing the auxiliary variable using the interpolated total mobility
                    intsj=intsj+log(rs/rf(u))*(lo/lt(u)-1+lo/ltskin-1)/2;
                    % if the next radius is inside the damaged zone, no interpolation is required
                else
                    % incrementing the auxiliary integral
                    intsj=intsj+log(rf(u-1)/rf(u))*(lo/lt(u)-1+lo/lt(u-1)-1)/2;
                end
                % decreasing the waterfront radius index until the skin radius is reached
                u = u - 1;
            end
            % multiplying the second integral by (kj/kskinj-1)
            intsj=intsj*(k/ks-1);
            % adding the second integral to aj
            aj=aj+intsj;
            % multiplying aj by alphap/(k*h*lohat)
            aj=aj/k/h/lo;
        end
    else
        % if there is no formation damage, multiplying aj by alphap/(k*h*lohat)
        aj=aj/k/h/lo;
    end
end

%% funcao que calcula o termo deltapo segundo o modelo de Lefkovits =======
function [x]=lefk(u,rw,q,h,k,eta,lo)
    % definindo o numero de camadas
    nlay=length(h);
    
    % inicializando x com zero
    x=0;
    % calculando a contribuicao de cada camada segundo eq 14 @pdf 5 do arquivo "02B-142746"
    for j=1:nlay
        % numerador da soma
        aux=lo(j)*k(j)*h(j)*rw*sqrt(u/eta(j))*besselk(1,sqrt(u/eta(j))*rw);
        % denominador da soma
        aux=aux/besselk(0,sqrt(u/eta(j))*rw);
        % incremetando a soma
        x=x+aux;
    end
    % invertendo a soma
    x=1/x;
    % multiplicando a soma pelo termo q/u
    x=x*q/u;
end

%% funcao que calcula o coef Vj, necessario para o algoritmo de Stehfest ==
function [vj]=Vj(N,j)
    % inicializando vj com zero
    vj=0.0;
    % calculando o valor de k1
    k1=floor((j+1)/2.0);
    % calculando o valor de kn
    kn=min(j,N/2);
    % implementando a soma
    for k=k1:kn
        % calculando o 1o fatorial
        a1=power(k,N*0.5)*factorial(2*k);
        % calculando o 2o fatorial
        a2=factorial(N*0.5-k);
        % calculando o 3o fatorial
        a3=factorial(k)*factorial(k-1);
        % calculando o 4o fatorial
        a4=factorial(j-k)*factorial(2*k-j);
        % incrementando a soma
        vj=vj+double(a1/a2/a3/a4);
    end
    % ajustando o sinal de vj
    vj=vj*power(-1,round(j+N/2));
end

%% funcao que calcula o raio da frente de avanco usando buckley leverett ==
function [rfxy] = calc_rf(ii, rw, t, q, h, phi, dfw)
    q=q/19.03;
    % computing the integral at the first timestep
    integral=q(1,:).*t(1);
    % numerically integrating (using the trapeze rule) the expression required to determine waterfront radius
    for u=1:ii-1
        integral=integral+(t(u+1)-t(u)).*(q(u+1,:)+q(u,:))./2;
    end
    % calculating the waterfront radii using eq. 54 from file "09-SPE-84957"
    rfxy=sqrt(integral.*dfw./24./pi./phi./h+rw*rw*1.01);
end

%% funcao que retorna os dados de permeabilidade relativa =================
function [sw,krw,kro] = get_krel(flag)
    switch flag
        case 0 % curva fornecida pelo Abelardo
            % inicializando a saturacao inicial de agua Swi
            swi=0.25;
            % inicializando a saturacao de oleo residual Sor
            sor=0.28;
            % definindo o numero de entradas na tabela de permeabilidade relativa
            compsw=21;
            % calculando o deltasw a partir do Swi, Sor e compsw
            deltasw=(1-sor-swi)/(compsw-1);

            % inicializando os vetores de saturacao e de permeabilidade relativa
            sw=zeros(compsw,1);
            krw=sw;
            kro=sw;
            
            % preenchendo o vetor de saturacoes (aqui, ele fica uniformemente espacado)
            sw=swi:deltasw:1-sor;
%             for ii=2:compsw
%                 sw(ii)=sw(ii-1)+deltasw;
%             end

            % dados de permeabilidade relativa da agua
            krw(1)=0.000000; krw(2)=0.004730; krw(3)=0.009458;
            krw(4)=0.014188; krw(5)=0.018916; krw(6)=0.023646;
            krw(7)=0.028374; krw(8)=0.033110; krw(9)=0.037873;
            krw(10)=0.042757; krw(11)=0.047613; krw(12)=0.052476;
            krw(13)=0.058981; krw(14)=0.066197; krw(15)=0.073700;
            krw(16)=0.082294; krw(17)=0.095461; krw(18)=0.110594;
            krw(19)=0.127831; krw(20)=0.148597; krw(21)=0.173000;

            % dados de permeabilidade relativa do oleo
            kro(1)=0.540300; kro(2)=0.500145; kro(3)=0.459990;
            kro(4)=0.419835; kro(5)=0.379685; kro(6)=0.339530;
            kro(7)=0.299375; kro(8)=0.261046; kro(9)=0.223727;
            kro(10)=0.188473; kro(11)=0.154791; kro(12)=0.122027;
            kro(13)=0.095422; kro(14)=0.071638; kro(15)=0.051815;
            kro(16)=0.034833; kro(17)=0.024335; kro(18)=0.015993;
            kro(19)=0.009207; kro(20)=0.003933; kro(21)=0.000000;

        case 1 % curva retirada de um video no youtube
                % inicializando a saturacao inicial de agua Swi
                swi=0.20;
                % inicializando a saturacao de oleo residual Sor
                sor=0.40;
                % definindo o numero de entradas na tabela de permeabilidade relativa
                compsw=33;
                % calculando o deltasw a partir do Swi, Sor e compsw
                deltasw=(1-sor-swi)/(compsw-1);
                % preenchendo o vetor de saturacoes (aqui, ele fica uniformemente espacado)
                sw=swi:deltasw:1-sor;

                % dados de permeabilidade relativa da agua
                krw(1)=0.00000000; krw(2)=0.0005208997; krw(3)=0.001857067;
                krw(4)=0.003906338; krw(5)=0.006620655; krw(6)=0.009968479;
                krw(7)=0.013926540; krw(8)=0.018476520; krw(9)=0.023603390;
                krw(10)=0.02929446; krw(11)=0.03553877; krw(12)=0.04232670;
                krw(13)=0.04964970; krw(14)=0.05750009; krw(15)=0.06587090;
                krw(16)=0.07475575; krw(17)=0.08414878; krw(18)=0.09404457;
                krw(19)=0.10443800; krw(20)=0.11532500; krw(21)=0.12670000;
                krw(22)=0.13855900; krw(23)=0.15090000; krw(24)=0.16371700;
                krw(25)=0.17700700; krw(26)=0.19076700; krw(27)=0.20499400;
                krw(28)=0.21968500; krw(29)=0.23483700; krw(30)=0.25044700;
                krw(31)=0.26651300; krw(32)=0.28303100; krw(33)=0.30000000;

                % dados de permeabilidade relativa do oleo
                kro(1)=0.800000000; kro(2)=0.754998000; kro(3)=0.711177000;
                kro(4)=0.668542000; kro(5)=0.627101000; kro(6)=0.586861000;
                kro(7)=0.547830000; kro(8)=0.510017000; kro(9)=0.473429000;
                kro(10)=0.438075000; kro(11)=0.403966000; kro(12)=0.371110000;
                kro(13)=0.339517000; kro(14)=0.309200000; kro(15)=0.280169000;
                kro(16)=0.252436000; kro(17)=0.226016000; kro(18)=0.200922000;
                kro(19)=0.177168000; kro(20)=0.154773000; kro(21)=0.133753000;
                kro(22)=0.114128000; kro(23)=0.095920380; kro(24)=0.079153240;
                kro(25)=0.063854010; kro(26)=0.050053620; kro(27)=0.037787920;
                kro(28)=0.027099420; kro(29)=0.018040030; kro(30)=0.010675840;
                kro(31)=0.005096669; kro(32)=0.001439911; kro(33)=0.000000000;

        case 2 % curva video modificada
            % dados de saturacao
            sw(1)=0.2000; sw(2)=0.2125; sw(3)=0.2375;
            sw(4)=0.2500; sw(5)=0.2750; sw(6)=0.2875;
            sw(7)=0.3250; sw(8)=0.3500; sw(9)=0.3625;
            sw(10)=0.3875; sw(11)=0.4000; sw(12)=0.4250;
            sw(13)=0.4375; sw(14)=0.4625; sw(15)=0.4750;
            sw(16)=0.5125; sw(17)=0.5375; sw(18)=0.5500;
            sw(19)=0.5750; sw(20)=0.5875; sw(21)=0.6000;

            % dados de permeabilidade relativa da agua
            krw(1)=0.00000; krw(2)=0.00052; krw(3)=0.00391;
            krw(4)=0.00662; krw(5)=0.01393; krw(6)=0.01848;
            krw(7)=0.03554; krw(8)=0.04965; krw(9)=0.05750;
            krw(10)=0.07476; krw(11)=0.08415; krw(12)=0.10444;
            krw(13)=0.11533; krw(14)=0.13856; krw(15)=0.15090;
            krw(16)=0.19077; krw(17)=0.21969; krw(18)=0.23484;
            krw(19)=0.26651; krw(20)=0.28303; krw(21)=0.30000;

            % dados de permeabilidade relativa do oleo
            kro(1)=0.80000; kro(2)=0.75500; kro(3)=0.66854;
            kro(4)=0.62710; kro(5)=0.54783; kro(6)=0.51002;
            kro(7)=0.40397; kro(8)=0.33952; kro(9)=0.30920;
            kro(10)=0.25244; kro(11)=0.22602; kro(12)=0.17717;
            kro(13)=0.15477; kro(14)=0.11413; kro(15)=0.09592;
            kro(16)=0.05005; kro(17)=0.02710; kro(18)=0.01804;
            kro(19)=0.00510; kro(20)=0.00144; kro(21)=0.00000;

        case 3 % deslocamento pistao
            % dados de saturacao
            sw(1)=0.0;  sw(2)=1.0;
            % dados de permeabilidade relativa da agua
            krw(1)=0.0; krw(2)=1.0;
            % dados de permeabilidade relativa do oleo
            kro(1)=1.0; kro(2)=0.0;

        case 4 % deslocamento pistao com mais pontos
            % definindo o numero de pontos
            compsw=10;
            
            % pontos extremais
            swi=0;
            sor=0;
            ko=1;
            kw=1;
            % deltas
            dsw=(1-swi-sor)/(compsw-1);
            dko=ko/(compsw-1);
            dkw=kw/(compsw-1);
            
            % dados de saturacao
            sw=swi:dsw:1-sor;
            % dados de kro
            kro=ko:-dko:0;
            % dados de krw
            krw=0:+dkw:kw;
            
        otherwise % deslocamento pistao
            % dados de saturacao
            sw(1)=0.0;  sw(2)=1.0;
            % dados de permeabilidade relativa da agua
            krw(1)=0.0; krw(2)=1.0;
            % dados de permeabilidade relativa do oleo
            kro(1)=1.0; kro(2)=0.0;
    end
end

%% funcao que recebe krel, muo e muw para calcular lt e dfw ==========
function [lt,dfw]=fill_data(compsw,krw,kro,sw,mio,miw)
    % inicializando com zeros os dados de saida e o vetor auxiliar fw
    fw=zeros(compsw,1);
    dfw=zeros(compsw,1);
    lt=zeros(compsw,1);
    
    % calculando os dados de mobilidade total
    lt(:)=kro./mio+krw./miw;
    % calculando os dados de fluxo fracionario
    fw(:)=(krw(:)./miw)./lt(:);
    
    % calculando numericamente a derivada do fluxo fracionario no 1o ponto
    dfw(1)=(fw(2)-fw(1))/(sw(2)-sw(1));
    % a derivada do fluxo fracionario no ultimo ponto deve ser sempre zero
    dfw(compsw)=0;
    % para cada ponto com 2 vizinhos, calculando a derivada usando Bourdet
    for ii=2:compsw-1
        % calculando os dois termos necessarios para a derivada de Bourdet
        a=((fw(ii+1,:)-fw(ii,:))/(sw(ii+1,:)-sw(ii,:)))*(sw(ii,:)-sw(ii-1,:));
        b=((fw(ii,:)-fw(ii-1,:))/(sw(ii,:)-sw(ii-1,:)))*(sw(ii+1,:)-sw(ii,:));
        % calculando a derivada de Bourdet
        dfw(ii,:)=(a+b)/(sw(ii+1,:)-sw(ii-1,:));
    end
    
    % usando o metodo de Welge para ajustar a derivada do fluxo fracionario
    % o metodo deve partir do final da curva
    index=compsw;
    % calculando a derivada auxiliar do fluxo fracionario
    dfwaux=(fw(index)-fw(1))/(sw(index)-sw(1));
    % enquanto a derivada auxiliar for maior que a derivada original,
    % continuar a busca pelo ponto de tangencia
    while (dfw(index)<dfwaux && index >1)
        % atualizando a derivada auxiliar do fluxo fracionario
        dfwaux=(fw(index)-fw(1))/(sw(index)-sw(1));
        % atualizando o indice no vetor de saturacoes
        index=index-1;
    end
    % atualizar o vetor original de derivada somente se o ponto de tangencia
    % nao for o primeiro ponto
    if (index>1)
        % substituindo todos os pontos de derivada do fluxo fracionario
        % abaixo do ponto de tangencia
        for ii=1:index
            dfw(ii)=dfwaux;
        end
    end
    
end

%% Funcao que calcula a derivada da pressao com relacao ao ln do tempo ====
function [dpwf]=compute_derivative(t,pwf,tp,flagperiodo)
    
    % definindo o tamanho do vetor
    dim=length(pwf);
    % initializando o vetor de derivada com zeros
    dpwf=zeros(size(pwf));
    
    % calculando a derivada no 1o ponto
    dpwf(1)=(pwf(2)-pwf(1))/log(t(2)/t(1));
    
    % encontrando o ponto de tempo em que t = tp
    if flagperiodo==2
        % se houver periodo de falloff, buscando o tempo em que t = tp
%         tflag=find(abs(t-tp)<1e-10);
        tflag=find(t>tp);
        tflag=tflag(1);
        % para evitar problemas de aproximacao numerica, a busca tolera um erro de 1e-12
        
        % para todos os pontos com pelo menos 2 vizinhos, calculando a
        % derivada segundo Bourdet
        for ii=2:tflag-1
            % calculando os 2 termos necessarios para obter a derivada de Bourdet
            a = (pwf(ii + 1) - pwf(ii)) / log(t(ii + 1) / t(ii))*log(t(ii) / t(ii - 1));
            b = (pwf(ii) - pwf(ii - 1)) / log(t(ii) / t(ii - 1))*log(t(ii + 1) / t(ii));
            % Derivada de Bourdet disponivel no pdf 74 do arquivo "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + 1) / t(ii-1));
        end
        % calculando a derivada no ultimo ponto
        dpwf(tflag)=(pwf(tflag)-pwf(tflag-1))/log(t(tflag)/t(tflag-1));
        
        % calculando os tempos equivalentes para o 1o ponto do falloff
        teq0=tp*(t(tflag+1)-tp)/t(tflag+1);
        teq1=tp*(t(tflag+2)-tp)/t(tflag+2);
        % calculando a derivada no primeiro ponto de falloff
        dpwf(tflag+1)=(pwf(tflag+2)-pwf(tflag+1))/log(teq1/teq0);
        % para os pontos com 2 vizinhos, calculando a derivada segundo Bourdet
        for ii=dim/2+2:dim-1
            % calculando os tempos equivalentes 
            teq0=tp*(t(ii-1)-tp)/t(ii-1);
            teq1=tp*(t(ii-0)-tp)/t(ii-0);
            teq2=tp*(t(ii+1)-tp)/t(ii+1);
            % calculando os 2 termos necessarios para obter a derivada de Bourdet
            a = (pwf(ii+1) - pwf(ii))*log(teq0 / teq1) / (log(teq2 / teq1)*log(teq2 / teq0));
            b = (pwf(ii) - pwf(ii-1))*log(teq1 / teq2) / (log(teq1 / teq0)*log(teq2 / teq0));
            % somando os 2 termos para obter a derivada de Bourdet
            dpwf(ii)=a+b;
        end
        % calculando os tempos equivalentes para o ultimo ponto do falloff
        teq1=tp*(t(dim-1)-tp)/t(dim-1);
        teq2=tp*(t(dim)-tp)/t(dim);
        % calculando a derivada no ultimo ponto do falloff
        dpwf(dim)=(pwf(dim)-pwf(dim-1))/log(teq1/teq2);
    else
        % se houver apenas o periodo de injecao, definindo tflag como a
        % ultima posicao do vetor de tempo
        tflag=length(t);
        
        % para todos os pontos com pelo menos 2 vizinhos, calculando a
        % derivada segundo Bourdet
        for ii=2:tflag-1
            % calculando os 2 termos necessarios para obter a derivada de Bourdet
            a = (pwf(ii + 1) - pwf(ii)) / log(t(ii + 1) / t(ii))*log(t(ii) / t(ii - 1));
            b = (pwf(ii) - pwf(ii - 1)) / log(t(ii) / t(ii - 1))*log(t(ii + 1) / t(ii));
            % Derivada de Bourdet disponivel no pdf 74 do arquivo "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + 1) / t(ii-1));
        end
        % calculando a derivada no ultimo ponto
        dpwf(tflag)=(pwf(tflag)-pwf(tflag-1))/log(t(tflag)/t(tflag-1));
    end
    
end
% =========================================================================
%   Created by Renan Vieira on 04/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   O codigo almeja calcular a pressao em um reservatorio com poco
%   vertical durante Tinj. Por enquanto, apenas 1 camada. Aqui,
%   tentamos aplicar a prod inst na CCI (e a EDP fica homog).
%   Nesse codigo, vamos tentar aplicar diferentes CCEs
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
fcce=3;

% definindo o tempo inicial
t0=9.6e-8;
% definindo o tempo de fechamento do poco (definido pelo usuario)
tp=9.6*1e1;
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
nlay=1;%length(fid(:,1));
% lendo o raio do poco
rw=fid(1,1);
% lendo a vazao
q=fid(1,2);

% lendo as permeabilidades na regiao 1
k1=fid(1,3)';
% lendo os raios da regiao 1
r1=fid(1,4)';
% lendo as permeabilidades na regiao 2
k2=fid(1,5)';
% lendo os raios da regiao 2
r2=fid(1,6)';
% lendo as permeabilidades na direcao z
kz2=fid(1,7)';
% lendo as porosidades das camadas
phi=fid(1,8)';
% lendo as espessuras das camadas
h=fid(1,9)';
% lendo o comprimento do poco em cada camada
len=fid(1,10)';
% lendo as viscosidades do oleo
muo=fid(1,11)';

%% definindo alguns parametros de entrada =================================
% definindo o valor da viscosidade da agua (em cP)
muw=0.5;
% definindo o fator volume de formacao da agua (em m^3/ STD m^3)
bw=1;
% buscando os flags de permeabilidade relativa (definidos pelo usuario na interface)
flag=fid(1,12)';

%% calculando parametros que sao definidos a partir dos dados de entrada ==

for jj=1:nlay
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(flag(jj));
end
dfw=zeros(length(kro(:,1)),nlay);
for jj=1:nlay
    [~,dfw(:,jj)]=fill_data(length(sw(:,jj)),krw(:,jj),kro(:,jj),sw(:,jj),muo(jj),muw);
end

lo=kro(1)./muo;
lw=krw(end)./muw;
% calculando a compressibilidade total (ct = cr + cw*swi + co*(1-sor))
ct=2.9e-5*ones(size(h)); % cr+sw(1).*cw+sw(end).*co; %
% calculando a difusividade hidraulica (eta = k2*lohat/(phi*ct))
etao=alphat.*k2*lo./phi./ct;
etaw=etao./lo.*lw;
% ajustando as unidades da vazao
q=alphap*q*bw;

%% sessao onde sera chamada a funcao que calcula a pressao e as vazoes ====
% calculando os dados de pressao e vazao
[pg, ~, rfg]=VW_gre6(t,nlay,q,rw,dfw,phi,len,k2,r2,etaw,etao,lw,lo,fcce);

% calculando os dados de pressao e vazao
[pa, ~, rfa]=VW_p(t,nlay,q,rw,dfw,phi,len,k2,r2,etaw,etao,lw,lo,3);

% calculando a derivada da pressao
[dpg]=compute_derivative(t,pg,tp,fper);
% calculando a derivada da pressao
[dpa]=compute_derivative(t,pa,tp,fper);

% calculando o patamar teorico da derivada 
mo=q*muo/2/k2/kro(1)/len;
mw=mo*kro(1)/muo*muw/krw(end);
% mw=qinj*miw/2/k/krw(end)/h;
clc;
%% plotando os resultados =================================================

% definindo os limites dos eixos nos graficos
bo=round(log10(t(1))); bo=max([10^(bo) 1e-4]); 
bw=floor(log10(tp))+1; bw=10^bw;

co=min(mo, mw); co=1*10^floor(log10(co));
cw=ceil(log10(max(pa))); cw=10^cw;
% co=cw/100;

linw=1.5;
fonts=14;
% plotando os dados durante o periodo de injecao
figure
% loglog(t,pwf,'-k','LineWidth',linw);
loglog(t,pg,'ok', t,dpg,'ob','LineWidth',linw);
hold on
loglog(t,pa,'-k', t,dpa,'-b','LineWidth',linw);
yline(mo,'--k')
yline(mw,':k','LineWidth',2*linw)
legend({'\Delta P (this work) ', 'der. (this work) ',...
    '\Delta P (Viana et al.) ', 'der. (Viana et al.) '},'Location',...
    'northwest', 'fontsize',fonts*.75)
% legend({'\Delta P (canal) ', 'der. (canal) ','\Delta P (no-flow) ', ...
%     'der. (no-flow) '},'Location','northwest', 'fontsize',fonts*.75)
axis([bo,bw, co,cw])
title('Pressure and Pressure Derivative Profile')
xlabel('t (h)', 'fontsize',fonts)
ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
grid on

q=q/alphap/bw;

%% limpando algumas variaveis que nao sao mais necessarias ================
clear alphap alphat bo bw co cr cw dim dfw fid filename flag fper jj kro
clear krw p0 linw fonts
% clear 
toc

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
function [pwf,qj,rf]=VW_gre6(t,nlay,q,rw,dfw,phi,h,k2,r2,eta1,eta2,l1,l2,fcce)
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay);
    
    % definindo o tamanho do infinito (usado para CCE canal)
    if fcce==3
        infi=90;
    end
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=2:dim
            qj(ii,jj)=q;
            % calculando o raio da frente de avanco
            rf=calc_rf(ii, rw, t, qj(:,jj)/1, h(jj), phi(jj), dfw(:,jj));
            rf=rf(1);
            
            % entrando no loop da convolucao
            for ia=2:ii
                imp1=0;
                imp2=0;
                imp3=0;
                % entrando no loop de stehfest
                for kk=1:N
                    % calculando a "variavel de Laplace"
                    u1=log(2)*kk/(t(ia-1));
                    u2=log(2)*kk/(t(ia+0)+t(ia-1))*2;
                    u3=log(2)*kk/(t(ia+0));
                    % calculando o valor de x
                    x1=Vj(N,kk)*sisgre(u1,1,h(jj),rw,rf,r2,l1,l2,eta1,eta2,k2,0,fcce);
                    x2=Vj(N,kk)*sisgre(u2,1,h(jj),rw,rf,r2,l1,l2,eta1,eta2,k2,0,fcce);
                    x3=Vj(N,kk)*sisgre(u3,1,h(jj),rw,rf,r2,l1,l2,eta1,eta2,k2,0,fcce);
                    % incrementando a soma
                    imp1=imp1+x1(1)*besseli(0,sqrt(rw*rw*u1/eta1));
                    imp1=imp1+x1(2)*besselk(0,sqrt(rw*rw*u1/eta1));
                    imp2=imp2+x2(1)*besseli(0,sqrt(rw*rw*u2/eta1));
                    imp2=imp2+x2(2)*besselk(0,sqrt(rw*rw*u2/eta1));
                    imp3=imp3+x3(1)*besseli(0,sqrt(rw*rw*u3/eta1));
                    imp3=imp3+x3(2)*besselk(0,sqrt(rw*rw*u3/eta1));
                    switch fcce
                        % somando os pocos imagem para canal
                        case 3
                            for n=1:infi
                                imp1=imp1+2*x1(3)*besselk(0,n*r2*sqrt(u1/eta2));
                                imp2=imp2+2*x2(3)*besselk(0,n*r2*sqrt(u2/eta2));
                                imp3=imp3+2*x3(3)*besselk(0,n*r2*sqrt(u3/eta2));
                            end
                            % somando o poco imagem para falha selante
                        case 4
                            imp1=imp1+x1(3)*besselk(0,r2*sqrt(u1/eta2));
                            imp2=imp2+x2(3)*besselk(0,r2*sqrt(u2/eta2));
                            imp3=imp3+x3(3)*besselk(0,r2*sqrt(u3/eta2));
                            % somando o poco imagem para aquifero lateral
                        case 5
                            imp1=imp1-x1(3)*besselk(0,r2*sqrt(u1/eta2));
                            imp2=imp2-x2(3)*besselk(0,r2*sqrt(u2/eta2));
                            imp3=imp3-x3(3)*besselk(0,r2*sqrt(u3/eta2));
                    end
                end
                if isnumeric(imp1)==false
                    imp1=0;
                end
                if isnumeric(imp2)==false
                    imp2=0;
                end
                if isnumeric(imp3)==false
                    imp3=0;
                end
                % invertendo a soma para o campo real
                imp1=imp1*log(2)/(t(ia-1));
                imp2=imp2*log(2)/(t(ia+0)+t(ia-1))*2;
                imp3=imp3*log(2)/(t(ia+0));
                % incrementando a integral
                pwf(ii)=pwf(ii)+(t(ia)-t(ia-1))*(imp1+4*imp2+imp3)/6;
            end
        end
    end
    
    pwf=pwf*q;
    
end

%% sistema linear do impulso instantaneo radial no espaco de Laplace ======
function [x]=sisgre(u,Q,len,rw,r1,r2,l1,l2,eta1,eta2,k,t0,fcce)
    
    % definindo os argumentos das funções de bessel
    arg1=sqrt(rw*rw*u/eta1);
    arg2=sqrt(r1*r1*u/eta1);
    arg3=sqrt(r1*r1*u/eta2);
    
    c=l2/l1*sqrt(eta1/eta2);
    
    % selecionando o sistema linear adequado com a CCE
    switch fcce
        case 1    % sistema linear para CPB
            arg4=sqrt(r2*r2*u/eta2);
            % inicializando o vetor com zeros
            b=zeros(4,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)
                0                0                  besseli(0,arg4)    besselk(0,arg4)];
        
        case 2    % sistema linear para NFB
            arg4=sqrt(r2*r2*u/eta2);
            % inicializando o vetor com zeros
            b=zeros(4,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)
                0                0                  besseli(1,arg4)   -besselk(1,arg4)];
        otherwise    % sistema linear para IAR (tambem usado em res canal)
            % inicializando o vetor com zeros
            b=zeros(3,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0
                besseli(0,arg2)   besselk(0,arg2)   -besselk(0,arg3)
                besseli(1,arg2)  -besselk(1,arg2)  c*besselk(1,arg3)];
    end
    % alterando a 1a entrada do vetor
    b(1)=-Q/len/l1/k/arg1*exp(-u*t0);
    
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
function [pwf, qj, rf]=VW_p(t,nlay,q,rw,dfw,phi,h,k2,r2,eta1,eta2,l1,l2,fcce)
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay);
    p2=pwf;
    
    % definindo o tamanho do infinito (usado para CCE canal)
    if fcce==3
        infi=100;
    end
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=1:dim
            qj(ii,jj)=q;
            % calculando o raio da frente de avanco
            rf=calc_rf(ii, rw, t, qj(:,jj), h(jj), phi(jj), dfw(:,jj));
            rf=rf(1);
            % entrando no loop de stehfest
            for kk=1:N
                % calculando a "variavel de Laplace"
                u=log(2)*kk/t(ii);
                % calculando o valor de x
                x=Vj(N,kk)*sis(u,q,h(jj),rw,rf,r2,l1,l2,eta1,eta2,k2,fcce);
                % incrementando a soma
                pwf(ii)=pwf(ii)+x(1)*besseli(0,sqrt(rw*rw*u/eta1))+...
                    x(2)*besselk(0,sqrt(rw*rw*u/eta1));
                switch fcce
                    % somando os pocos imagem para canal
                    case 3
                        for n=1:infi
                            pwf(ii)=pwf(ii)+2*x(3)*besselk(0,n*r2*sqrt(u/eta2));
                        end
                % somando o poco imagem para falha selante
                    case 4
                    pwf(ii)=pwf(ii)+x(3)*besselk(0,sqrt(r2*r2*u/eta2));
                % somando o poco imagem para aquifero lateral
                    case 5
                    pwf(ii)=pwf(ii)-x(3)*besselk(0,sqrt(r2*r2*u/eta2));
                end
            end
            % invertendo a soma para o campo real
            pwf(ii)=pwf(ii)*log(2)/t(ii);
            p2(ii)=p2(ii)*log(2)/t(ii);
        end
    end
    
end

%% sistema linear do modelo radial composto no espaco de Laplace ==========
function [x]=sis(u,q,len,rw,r1,r2,l1,l2,eta1,eta2,k,fcce)

    % definindo os argumentos das funções de bessel
    arg1=sqrt(rw*rw*u/eta1);
    arg2=sqrt(r1*r1*u/eta1);
    arg3=sqrt(r1*r1*u/eta2);

    c=l2/l1*sqrt(eta1/eta2);
    
    switch fcce
        case 1    % sistema linear para CPB
            arg4=sqrt(r2*r2*u/eta2);
            % inicializando o vetor com zeros
            b=zeros(4,1);
            b(4)=1*0/u;
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)
                0                0                  besseli(0,arg4)    besselk(0,arg4)];
        
        case 2    % sistema linear para NFB
            arg4=sqrt(r2*r2*u/eta2);
            % inicializando o vetor com zeros
            b=zeros(4,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0   0
               besseli(0,arg2)   besselk(0,arg2)   -besseli(0,arg3)   -besselk(0,arg3)
               besseli(1,arg2)  -besselk(1,arg2) -c*besseli(1,arg3)  c*besselk(1,arg3)
                0                0                  besseli(1,arg4)   -besselk(1,arg4)];
        
        otherwise    % sistema linear para IAR (tambem usado na CCE canal)
            % inicializando o vetor com zeros
            b=zeros(3,1);
            
            % calculando a matriz a ser invertida
            m=[besseli(1,arg1)  -besselk(1,arg1)	0
                besseli(0,arg2)   besselk(0,arg2)   -besselk(0,arg3)
                besseli(1,arg2)  -besselk(1,arg2)  c*besselk(1,arg3)];
    end
    % alterando a 1a entrada do vetor
    b(1)=-q/u/len/k/l1/arg1;
    
    a=max(max(m));
    m=m/a;
    b=b/a;

    % calculando os coeficientes
    x=m\b;

%     % calculando a pressao no poco
%     y=x(1)*besseli(0,arg1)+x(2)*besselk(0,arg1);
% 
%     % conferindo se a continuidade foi mantida
%     ya=x(1)*besseli(0,arg2)+x(2)*besselk(0,arg2);
%     yx=x(3)*besselk(0,arg3);
%     disp(ya-yx)
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
    rfxy=sqrt(integral.*dfw./24./pi./phi./h+rw*rw);
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
            for ii=2:compsw
                sw(ii)=sw(ii-1)+deltasw;
            end

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

        otherwise % deslocamento pistao
            % dados de saturacao
            sw(1)=0.0;  sw(2)=1.0;
            % dados de permeabilidade relativa da agua
            krw(1)=0.0; krw(2)=1.0;
            % dados de permeabilidade relativa do oleo
            kro(1)=1.0; kro(2)=0.0;
    end
end

%% funcao que recebe krel, muo e muw para calcular lambdat e dfw ==========
function [lambdat,dfw]=fill_data(compsw,krw,kro,sw,mio,miw)
    % inicializando com zeros os dados de saida e o vetor auxiliar fw
    fw=zeros(compsw,1);
    dfw=zeros(compsw,1);
    lambdat=zeros(compsw,1);
    
    % calculando os dados de mobilidade total
    lambdat(:)=kro./mio+krw./miw;
    % calculando os dados de fluxo fracionario
    fw(:)=(krw(:)./miw)./lambdat(:);
    
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
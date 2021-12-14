% =========================================================================
%   Created by Renan Vieira on 04/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   O codigo almeja calcular a pressao em um reservatorio com poco
%   vertical e onde ha regioes concentricas com o poco com 
%   mobilidades distintas. Por enquanto, monof e apenas 1 camada.
%   Aqui, tentamos aplicar a prod inst na CCI (e a EDP fica homog).
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

% definindo o tempo inicial
t0=9.6e-5;
% definindo o tempo de fechamento do poco (definido pelo usuario)
tp=96;
% definindo o numero de passos em cada periodo 
dim=round(log10(tp/t0));
dim=10*dim+1;

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
filename='props rad.txt';
fid=importdata(filename);
fid=fid.data;

% determinando o numero de camadas a partir do arquivo txt 
nlayers=1;%length(fid(:,1));
% lendo o raio do poco
rw=fid(1,10);
% lendo a vazao
qinj=fid(1,11);
% ajustando as unidades da vazao
qinj=alphap*qinj;

% lendo as porosidades das camadas
phi=fid(1,1)';
% lendo as permeabilidades na regiao 2 (por enquanto corresponde ao kx)
k2=fid(1,2)';
% lendo os raios da regiao 1 (por enquanto corresponde ao rskin)
r1=fid(1,4)';
% lendo as permeabilidades na regiao 1 (por enquanto corresponde ao kskin)
k1=fid(1,5)';
% lendo o comprimento do poco em cada camada
len=fid(1,6)';
% lendo as espessuras das camadas
h=fid(1,7)';
% lendo as viscosidades do oleo
mio=fid(1,8)';

%% definindo alguns parametros de entrada =================================
% definindo o valor da viscosidade da agua (em cP)
% miw=0.5;
% definindo o fator volume de formacao da agua (em m^3/ STD m^3)
bw=1;

%% calculando parametros que sao definidos a partir dos dados de entrada ==
l2=k2./mio;
l1=k1./mio;
% calculando a compressibilidade total (ct = cr + cw*swi + co*(1-sor))
ct=1.1e-4*ones(size(h)); % cr+sw(1).*cw+sw(end).*co; %
% calculando a difusividade hidraulica (eta = k2*lohat/(phi*ct))
eta2=alphat.*l2./phi./ct;
eta1=eta2./l2.*l1;

%% sessao onde sera chamada a funcao que calcula a pressao e as vazoes ====
% calculando os dados de pressao e vazao
[pg, ~]=HW_pres2(t,nlayers,qinj,rw,bw,h,eta1,eta2,l1,l2,r1);

% calculando os dados de pressao e vazao
[pa, ~]=HW_p(t,nlayers,qinj,rw,bw,h,eta1,eta2,l1,l2,r1);

% calculando a derivada da pressao
[dpg]=compute_derivative(t,pg,tp,fper);
% calculando a derivada da pressao
[dpa]=compute_derivative(t,pa,tp,fper);

% calculando o patamar teorico da derivada 
m2=qinj*mio/2/k2/h;
m1=m2*k2/k1;
% mw=qinj*miw/2/k2/krw(end)/h;

%% plotando os resultados =================================================

% definindo os limites dos eixos nos graficos
bo=round(log10(t(1))); bo=max([10^(bo) 1e-4]); 
bw=floor(log10(tp))+1; bw=10^bw;

co=min(m2, m1); co=10^floor(log10(co));
cw=floor(log10(max(pg)))+1; cw=10^cw;
% co=cw/100;

linw=1.5;
fonts=14;
% plotando os dados durante o periodo de injecao
figure
% loglog(t,pwf,'-k','LineWidth',linw);
loglog(t,pg,'vk', t,dpg,'vb','LineWidth',linw);
hold on
loglog(t,pa,'-k', t,dpa,'-b','LineWidth',linw);
yline(m2,'--k')
yline(m1,':k','LineWidth',linw)
% legend({'P Green ', 'dP Green ','P ant. ', 'dP ant. ', 'mo', 'mw'},...
%     'Location','northwest', 'fontsize',fonts*.85)
legend({'P Green ', 'dP Green ','P Isa ', 'dP Isa '},...
    'Location','northwest', 'fontsize',fonts*.85)
axis([bo,bw, co,cw])
title('Pressure and Pressure Derivative Profile')
xlabel('t (h)', 'fontsize',fonts)
ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
grid on

qinj=qinj/alphap;

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
function [pwf, qj]=HW_pres3(t,nlay,q,rw,bw,h,eta1,eta2,l1,l2,r1)
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay)+0*norm(h);
    
    % inicializando o vetor que armazena os impulsos radiais
%     impr=pwf;
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
        
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=2:dim
            qj(ii,jj)=q;
            % entrando no loop da convolucao
            for ia=2:ii
                imp1=0;
                imp2=0;
                imp3=0;
                % entrando no loop de stehfest
                for kk=1:N
                    % calculando a "variavel de Laplace"
                    u1=log(2)*kk/t(ia-1);
%                     u2=log(2)*kk/(t(ia+0)-t(ia-1))*2;
                    u3=log(2)*kk/t(ia+0);
                    % calculando o valor de x
                    x1=Vj(N,kk)*sis2(u1,1,h(jj),rw,r1,l1,l2,eta1,eta2,0);
%                     x2=Vj(N,kk)*sis2(u2,1,h(jj),rw,r1,l1,l2,eta1,eta2,0);
                    x3=Vj(N,kk)*sis2(u3,1,h(jj),rw,r1,l1,l2,eta1,eta2,0);
                    % incrementando a soma
                    imp1=imp1+x1;
%                     imp2=imp2+x2;
                    imp3=imp3+x3;
                end
                % invertendo a soma para o campo real
                imp1=imp1*log(2)/t(ia-1);
                imp2=imp2*log(2)/(t(ia+0)+t(ia-1))*2;
                imp3=imp3*log(2)/t(ia+0);
                % incrementando a integral
%                 pwf(ii)=pwf(ii)+(t(ia)-t(ia-1))*(imp1+4*imp2+imp3)/6;
                pwf(ii)=pwf(ii)+(t(ia)-t(ia-1))*(imp1+4*imp2+imp3)/2;
            end
        end
    end
    
    pwf=pwf*q*bw;
    
end

%% nova funcao que calcula a pressao a cada passo de tempo ================
function [pwf, qj]=HW_pres2(t,nlay,q,rw,bw,h,eta1,eta2,l1,l2,r1)
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay)+0*norm(h);
    
    % inicializando o vetor que armazena os impulsos radiais
    impr=pwf;
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % determinando o volume produzido em cada producao instantanea
    c=1/24;
    Q=q*bw*t(dim)*c;
    Q=Q/dim;
    Q=Q/Q;
    
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=1:dim
            qj(ii,jj)=q;
            % entrando no loop de stehfest
            for kk=1:N
                % calculando a "variavel de Laplace"
                u=log(2)*kk/t(ii);
%                 u=log(2)*kk/t(dim-ii+1);
                % calculando o valor de x
                x=Vj(N,kk)*sis2(u,Q,h(jj),rw,r1,l1,l2,eta1,eta2,0);
                % incrementando a soma
                impr(ii)=impr(ii)+x;
            end
            % invertendo a soma para o campo real
            impr(ii)=impr(ii)*log(2)/t(ii);
        end
    end
    
    % integrando as producoes pontuais
    pwf(1)=impr(1)*t(1)*4;
    for ii=2:dim
%         pwf(ii)=pwf(ii-1)+pwf(ii)*(t(ii)-t(ii-1));
        pwf(ii)=pwf(ii-1)+(impr(ii-1)+impr(ii))*(t(ii)-t(ii-1))/2;
    end
    
    pwf=pwf*q*bw;
    
end

%% novo sistema linear do impulso instantaneo radial no espaco de Laplace =
function [y]=sis2(u,Q,len,rw,r1,l1,l2,eta1,eta2,t0)
    
    % definindo os argumentos das funções de bessel
    arg1=sqrt(rw*rw*u/eta1);
    arg2=sqrt(r1*r1*u/eta1);
    arg3=sqrt(r1*r1*u/eta2);
    
    c=l2/l1*sqrt(eta1/eta2);
    
    % inicializando o vetor com zeros
    b=zeros(3,1);
    % alterando a 1a entrada do vetor
    b(1)=-Q/len/l1/arg1*exp(-u*t0);
    
    % calculando a matriz a ser invertida
    m=[besseli(1,arg1)  -besselk(1,arg1)	0
       besseli(0,arg2)   besselk(0,arg2)   -besselk(0,arg3)
       besseli(1,arg2)  -besselk(1,arg2)  c*besselk(1,arg3)];
   
   a=max(max(m));
   m=m/a/1e20;
   b=b/a/1e20;
   % calculando os coeficientes
%    x=b\m;
%    x=inv(m)*b;
   x=m\b;
   
   % calculando a pressao no poco
   y=x(1)*besseli(0,arg1)+x(2)*besselk(0,arg1);
   
   % conferindo se a continuidade foi mantida
%    ya=x(1)*besseli(0,arg2)+x(2)*besselk(0,arg2);
%    yx=x(3)*besselk(0,arg3);
%    disp(ya-yx)
end

%% funcao que calcula a pressao pelo modelo rad comp ======================
function [pwf, qj]=HW_p(t,nlay,q,rw,bw,h,eta1,eta2,l1,l2,r1)
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay)+0*norm(h);
    
    q=q*bw;
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=1:dim
            qj(ii,jj)=q;
            % entrando no loop de stehfest
            for kk=1:N
                % calculando a "variavel de Laplace"
                u=log(2)*kk/t(ii);
                % calculando o valor de x
                x=Vj(N,kk)*sis_comp(u,q,h(jj),rw,r1,l1,l2,eta1,eta2);
                % incrementando a soma
                pwf(ii)=pwf(ii)+x;
            end
            % invertendo a soma para o campo real
            pwf(ii)=pwf(ii)*log(2)/t(ii);
        end
    end
    
end

%% sistema linear do modelo radial composto no espaco de Laplace ==========
function [y]=sis_comp(u,q,len,rw,r1,l1,l2,eta1,eta2)
    
    % definindo os argumentos das funções de bessel
    arg1=sqrt(rw*rw*u/eta1);
    arg2=sqrt(r1*r1*u/eta1);
    arg3=sqrt(r1*r1*u/eta2);
    
    c=l2/l1*sqrt(eta1/eta2);
    
    % inicializando o vetor com zeros
    b=zeros(3,1);
    % alterando a 1a entrada do vetor
    b(1)=-q/u/len/l1/arg1;
    
    % calculando a matriz a ser invertida
    m=[besseli(1,arg1)  -besselk(1,arg1)	0
       besseli(0,arg2)   besselk(0,arg2)   -besselk(0,arg3)
       besseli(1,arg2)  -besselk(1,arg2)  c*besselk(1,arg3)];
   
   a=max(max(m));
   m=m/a;
   b=b/a;
   
   % calculando os coeficientes
   x=m\b;
   
   % calculando a pressao no poco
   y=x(1)*besseli(0,arg1)+x(2)*besselk(0,arg1);
   
   % conferindo se a continuidade foi mantida
   ya=x(1)*besseli(0,arg2)+x(2)*besselk(0,arg2);
   yx=x(3)*besselk(0,arg3);
   disp(ya-yx)
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

% =========================================================================
%   Created by Renan Vieira on 04/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   O codigo almeja calcular a pressao em um reservatorio com poco
%   vertical durante well test. Por enquanto, apenas 1 camada. 
%   Aqui, tentamos aplicar a prod inst na CCI (e a EDP fica homog).
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

% flag da CCE (0 = IAR; 1 = CPB; 2 = NFB; 3 = Cha; 4 = SeF; 5 = LAq)
fcce=3;

% definindo o tempo inicial
t0=9.6e-9;
% definindo o tempo de fechamento do poco (definido pelo usuario)
tp=9.6*1e+2;
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
% definindo arbitrariamente o re
re=2*r2;
% lendo as permeabilidades regiao 3 (aproveitei o kz3)
k3=fid(1,7)';
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
% miw=0.5;
% definindo o fator volume de formacao da agua (em m^3/ STD m^3)
bw=1;

%% calculando parametros que sao definidos a partir dos dados de entrada ==
l3=k3./muo;
l2=k2./muo;
l1=k1./muo;
% calculando a compressibilidade total (ct = cr + cw*swi + co*(1-sor))
ct=2.9e-5*ones(size(h)); % cr+sw(1).*cw+sw(end).*co; %
% calculando a difusividade hidraulica (eta = k2*lohat/(phi*ct))
eta3=alphat.*l3./phi./ct;
eta2=alphat.*l2./phi./ct;
eta1=eta2./l2.*l1;
% ajustando as unidades da vazao
q=alphap*q*bw;

%% sessao onde sera chamada a funcao que calcula a pressao e as vazoes ====
% calculando os dados de pressao e vazao
[pg, ~]=VW_pgre(t,nlay,q,rw,h,eta1,eta2,eta3,l1,l2,l3,r1,r2,re,fcce);

% calculando os dados de pressao e vazao
[pa, ~]=VW_p(t,nlay,q,rw,h,eta1,eta2,eta3,l1,l2,l3,r1,r2,re,fcce);

% calculando a derivada da pressao
[dpg]=compute_derivative(t,pg,tp,fper);
% calculando a derivada da pressao
[dpa]=compute_derivative(t,pa,tp,fper);

% calculando o patamar teorico da derivada 
m2=q*muo/2/k2/h;
m1=m2*k2/k1;
% mw=qinj*miw/2/k2/krw(end)/h;

%% plotando os resultados =================================================

% definindo os limites dos eixos nos graficos
bo=round(log10(t(1))); bo=max([10^(bo) 1e-4]); 
bw=floor(log10(tp))+1; bw=10^bw;

co=min(m2, m1); co=10^floor(log10(co));
cw=ceil(log10(max(pa))); cw=10^cw;

linw=1.5;
fonts=14;
% plotando os dados durante o periodo de injecao
figure
% loglog(t,pwf,'-k','LineWidth',linw);
loglog(t,pg,'ok', t,dpg,'ob','LineWidth',linw);
hold on
loglog(t,pa,'-k', t,dpa,'-b','LineWidth',linw);
yline(m2,'--k')
yline(m1,':k','LineWidth',linw)
legend({'\Delta P (this work) ', 'der. (this work) ',...
    '\Delta P (Viana et al.) ', 'der. (Viana et al.) '},'Location',...
    'northwest', 'fontsize',fonts*.75)
axis([bo,bw, co,cw])
title('Pressure and Pressure Derivative Profile')
xlabel('t (h)', 'fontsize',fonts)
ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
grid on

q=q/alphap;

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
function [pwf, qj]=VW_pgre(t,nlay,q,rw,h,eta1,eta2,eta3,l1,l2,l3,r1,r2,re,fcce)
    iaux=find(t>1e-5); iaux=iaux(1);
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay)+0*norm(h);
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % definindo o tamanho do infinito (usado para CCE canal)
    if fcce==3
        infi=100;
    end
    
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=iaux:dim
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
                    u2=log(2)*kk/(t(ia+0)-t(ia-1))*2;
                    u3=log(2)*kk/t(ia+0);
                    % calculando o valor de x
                    x1=Vj(N,kk)*sisgre(u1,1,h(jj),rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,0,fcce);
                    x2=Vj(N,kk)*sisgre(u2,1,h(jj),rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,0,fcce);
                    x3=Vj(N,kk)*sisgre(u3,1,h(jj),rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,0,fcce);
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
                                imp1=imp1+2*x1(5)*besselk(0,n*r2*sqrt(u1/eta2));
                                imp2=imp2+2*x2(5)*besselk(0,n*r2*sqrt(u2/eta2));
                                imp3=imp3+2*x3(5)*besselk(0,n*r2*sqrt(u3/eta2));
                            end
                            % somando o poco imagem para falha selante
                        case 4
                            imp1=imp1+x1(5)*besselk(0,r2*sqrt(u1/eta2));
                            imp2=imp2+x2(5)*besselk(0,r2*sqrt(u2/eta2));
                            imp3=imp3+x3(5)*besselk(0,r2*sqrt(u3/eta2));
                            % somando o poco imagem para aquifero lateral
                        case 5
                            imp1=imp1-x1(5)*besselk(0,r2*sqrt(u1/eta2));
                            imp2=imp2-x2(5)*besselk(0,r2*sqrt(u2/eta2));
                            imp3=imp3-x3(5)*besselk(0,r2*sqrt(u3/eta2));
                    end
                end
                if isnan(imp1)==true
                    imp1=0;
                end
                if isnan(imp2)==true
                    imp2=0;
                end
                if isnan(imp3)==true
                    imp3=0;
                end
                % invertendo a soma para o campo real
                imp1=imp1*log(2)/t(ia-1);
                imp2=imp2*log(2)/(t(ia+0)+t(ia-1))*2;
                imp3=imp3*log(2)/t(ia+0);
                % incrementando a integral
                pwf(ii)=pwf(ii)+(t(ia)-t(ia-1))*(imp1+4*imp2+imp3)/6;
%                 pwf(ii)=pwf(ii)+(t(ia)-t(ia-1))*(imp1+4*imp2+imp3)/2;
            end
        end
    end
    
    pwf=pwf*q;
    
end

%% novo sistema linear do impulso instantaneo radial no espaco de Laplace =
function [x]=sisgre(u,Q,len,rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,t0,fcce)
    
    % definindo os argumentos das funções de bessel
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
    b(1)=-Q/len/l1/arg1*exp(-u*t0);
    
    a=max(max(m));
    m=m/a;
    b=b/a;
    % calculando os coeficientes
    %    x=b\m;
    %    x=inv(m)*b;
    x=m\b;
    
    % calculando a pressao no poco
%     y=x(1)*besseli(0,arg1)+x(2)*besselk(0,arg1);
    
end

%% funcao que calcula a pressao pelo modelo rad comp ======================
function [pwf, qj]=VW_p(t,nlay,q,rw,h,eta1,eta2,eta3,l1,l2,l3,r1,r2,re,fcce)
    
    % obtendo a dimensao do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressao e a matriz de vazoes com zeros
    pwf = zeros(dim,1);
    qj = zeros(dim,nlay)+0*norm(h);
    
    % definindo o tamanho do infinito (usado para CCE canal)
    if fcce==3
        infi=100;
    end
    
    % definindo o numero N para o algoritmo de Stehfest
    N=12;
    
    % para cada camada, calcular as respectivas fontes instantaneas
    for jj=1:nlay
        % para cada tempo, calcular as respectivas fontes instantaneas
        for ii=2:dim
            qj(ii,jj)=q;
            % entrando no loop de stehfest
            for kk=1:N
                % calculando a "variavel de Laplace"
                u=log(2)*kk/t(ii);
                % calculando o valor de x
                x=Vj(N,kk)*sis(u,q,h(jj),rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,fcce);
                % incrementando a soma
                pwf(ii)=pwf(ii)+x(1)*besseli(0,sqrt(rw*rw*u/eta1))+...
                    x(2)*besselk(0,sqrt(rw*rw*u/eta1));
                switch fcce
                    % somando os pocos imagem para canal
                    case 3
                        for n=1:infi
                            pwf(ii)=pwf(ii)+2*x(5)*besselk(0,n*r2*sqrt(u/eta2));
                        end
                % somando o poco imagem para falha selante
                    case 4
                    pwf(ii)=pwf(ii)+x(5)*besselk(0,sqrt(r2*r2*u/eta2));
                % somando o poco imagem para aquifero lateral
                    case 5
                    pwf(ii)=pwf(ii)-x(5)*besselk(0,sqrt(r2*r2*u/eta2));
                end
            end
            % invertendo a soma para o campo real
            pwf(ii)=pwf(ii)*log(2)/t(ii);
        end
    end
    
end

%% sistema linear do modelo radial composto no espaco de Laplace ==========
function [x]=sis(u,q,len,rw,r1,r2,re,l1,l2,l3,eta1,eta2,eta3,fcce)

    % definindo os argumentos das funções de bessel
    arg1=sqrt(rw*rw*u/eta1);
    arg2=sqrt(r1*r1*u/eta1);
    arg3=sqrt(r1*r1*u/eta2);
    arg4=sqrt(r2*r2*u/eta2);
    arg5=sqrt(r2*r2*u/eta3);
    
    c=l2/l1*sqrt(eta1/eta2);
    d=l3/l2*sqrt(eta2/eta3);
    
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
    b(1)=-q/u/len/l1/arg1;
    
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

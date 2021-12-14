% =========================================================================
%   Criado por Renan Vieira em 04/21.
%   O script aplica o algortimo de Nelder Mead para estimar as 
%   permeabilidades de cada camada considerando producao por um 
%   poco vertical em reservatorio multicamada homogeneo sem skin e
%   sem estocagem. Pode-se testar 2 funcoes objetivo: uma considera
%   somente dados de pressao e outra considera pressao e vazoes.
% =========================================================================

% limpando o console, a memoria e fechando as janelas graficas
clc; clear; close all

%% definindo alguns parametros de entrada
% definindo o flag da funcao objetivo (1 = somente dados de pressao,
% 2 = dados de pressao e vazao)
flag=2;

% definindo o numero de camadas (dimensao do problema)
nlay=3;
% definindo o limite inferior em x1
x1=1*ones(nlay,1);
% definindo o limite superior em x1
x2=3*ones(nlay,1);

% definindo a tolerancia
tol=1e-10;
% definindo o numero maximo de iteracoes
nmax=1500;

%% definindo as propriedades do reservatorio

% viscosidade do oleo
muo=5.1;
% vazao de producao
q=100;
% raio do poco
rw=0.1;
% espessura das camadas
h=15*ones(1,nlay);
% porosidade (igual em todas as camadas)
phi=0.2;
% compressibilidade total (igual em todas as camadas)
ct=1e-5;

% definindo as permeabilidades reais
kref=rand(1,nlay)*(x2(1)-x1(1))+x1(1);
kref=10.^(kref);

%% calculando o vetor de tempo
% tempo inicial
t0=1e-4;
% razao da PG usada para criar o vetor de tempo
dt=10^(0.1);
% dimensao do vetor de tempo
dim=61;
% inicializando o vetor de tempo com zeros
t=zeros(dim,1);
% preenchendo o vetor de tempo seguindo os termos de uma PG
t(1)=t0;
for ii=2:dim
    t(ii)=t(ii-1)*dt;
end

%% calculando as pressoes e vazoes de referencia
[ptrue, qjtrue]=calc_pres(t, nlay, phi, kref, ct, h, muo, rw, q);

%% adicionando um ruido gaussiano no dado de referencia
% gerando vetores aleatorios usando uma gaussiana de media zero e sdv 1
erp=randn(dim,1);
erq=randn(dim,nlay);

% definindo o percentual de ruido a ser aplicado
err=0.05;
% definindo o dado obserdado como sendo o dado original mais um ruido
pobs=ptrue+err*min(ptrue)*erp;
qjobs=zeros(size(qjtrue));
for j=1:nlay
    qjobs(:,j)=qjtrue(:,j)+err*min(qjtrue(:,j))*erq(:,j);
end

%% calculando o minimo usando o algoritmo de Nelder-Mead

[x, xk, tNM, res]=NM_search(nlay,x1,x2,tol,nmax,flag,t,phi,ct,h,muo,rw,q,pobs,qjobs); 

% desfazendo a transformacao logaritmica
kest=10.^x';

% exibindo os resultados
fprintf(' <strong> Algoritmo de Nelder-Mead </strong> \n')
fprintf('Tempo decorrido: %.4f      No iteracoes: %d\n',tNM,length(res))
fprintf('Chute inicial:\n')
kaux=10.^xk(:,1);
for ii=1:nlay
    fprintf('k%d: %.1f	',ii,kaux(ii))
end
[paux,qjaux]=calc_pres(t,nlay,phi,kaux',ct,h,muo,rw,q);
aux=fobj(paux,qjaux,pobs,qjobs,flag);
fprintf('\n Valor da funcao em Xini: %.2f\n',aux)
fprintf('Ponto estimado:\n')
for ii=1:nlay
    fprintf('k%d: %.1f	',ii,kest(ii))
end
[pest,qjest]=calc_pres(t,nlay,phi,kest,ct,h,muo,rw,q);
aux=fobj(pest,qjest,pobs,qjobs,flag);
fprintf('\n Valor da funcao em Xest: %.3f\n',aux)

fprintf(' <strong> Permeabilidades de Referencia </strong> \n')
for ii=1:nlay
    fprintf('k%d: %.1f	',ii,kref(ii))
end
fprintf('\n')

% calculando o valor da funcao a cada iteracao
zest=zeros(1,size(xk,2));
for jj=1:size(xk,2)
    kaux=10.^xk(:,jj);
    [paux,qjaux]=calc_pres(t,nlay,phi,kaux',ct,h,muo,rw,q);
    zest(jj)=fobj(paux,qjaux,pobs,qjobs,flag);
end

% plotando o residuo em cada iteracao
figure
yyaxis left
semilogy(res,'Linewidth',1.0)
yline(tol,'k--')
ylabel('Valor do residuo')
xlabel('iteracao')
yyaxis right
semilogy(zest,'Linewidth',1.0)
ylabel('Funcao avaliada em xk')

figure
semilogx(t,qjest,t,qjobs,'v','LineWidth',1.5)
ylim([0 1.1*max(max(qjobs))])
ylabel('flow-rate (m³/d)')
xlabel('t (h)')
figure
loglog(t,pest,'-b',t,pobs,'vr','LineWidth',1.5)
ylabel('\Delta P (kgf/cm²)')
xlabel('t (h)')

%% limpando algumas variaveis auxiliares
clear ii jj delta zref mx my s alphap alphat nt t0 dt
clear dim erp erq j x1 x2 zest paux qjaux aux kaux ptrue qjtrue x

%% algoritmo de Nelder-Mead
function [xbar, xk, t, res]=NM_search(n,xlow,xhigh,tol,imax,flag,t,phi,ct,h,muo,rw,q,pref,qref)
    % iniciando a contagem de tempo que leva para o algoritmo convergir
    tic
    
    % definindo o coeficiente de reflexao
    ar=1.0;
    % definindo o coeficiente de expansao
    ae=2.0;
    % definindo o coeficiente de contracao
    ac=0.5;
    % definindo o coeficiente de encolhimento
    as=0.5;    
    
    % inicializando com zeros a matriz que armazena xbar a cada iteracao
    xk=zeros(n,imax);
    % inicializando o vetor de residuos com zeros
    res=zeros(1,imax);
    % inicializando o simplex inicial com zeros
    x=zeros(n,n+1);
    % inicializando o vetor que armazena o valor da funcao em cada vertice
    y=zeros(1,n+1);
    
    % definindo o simplex inicial com valores aleatorios entre xlow e xhigh
    for ii=1:n
        x(ii,:)=rand(1,n+1)*(xhigh(ii)-xlow(ii))+xlow(ii);
    end
    
    % inicializando o contador de iteracoes
    it=1;
    % aqui foi usado outro criterio de convergencia, que depende da funcao
    % objetivo na iteracao anterior e na iteracao atual. Entao nessa linha
    % a funcao obj na iteracao anterior foi inicializada com um valor alto 
    ypre=100;
    
    % entrando no loop do algoritmo
    while true
        k=10.^x;
        % avaliando o valor da funcao em cada vertice do simplex inicial
        for jj=1:n+1
            [p,qj]=calc_pres(t,n,phi,k(:,jj)',ct,h,muo,rw,q);
            y(jj)=fobj(p,qj,pref,qref,flag);
        end
        
        % encontrando o indice correspondente ao menor valor de y
        jl=find(y==min(y));
        jl=jl(1);
        % armazenando o minimo da funcao para o simplex atual
        yl=y(jl);
        % encontrando o indice correspondente ao maior valor de y
        jh=find(y==max(y));
        jh=jh(1);
        % encontrando o vertice que apresenta o menor valor
        xh=x(:,jh);
        % armazenando o maximo da funcao para o simplex atual
        yh=y(jh);
        
        % inicializando o centroide com zeros
        xbar=zeros(n,1);
        % calculando o centroide de todos os vertices exceto xh
        for jj=1:n+1
            % somando todos os vertices
            xbar=xbar+x(:,jj);
        end
        % calculando a media sem xh
        xbar=(xbar-xh)/n;
        
        % atualizando o residuo (criterio de convergencia 0)
        res(it)=abs(ypre-yh);
        
        % incluindo a estimativa atual na matriz de iteracoes
        xk(:,it)=xbar;
        
        % verificando se o numero maximo de iteracoes foi excedido
        % ou se o residuo esta menor que a tolerancia estabelecida
        if (it >= imax || res(it) < tol)
            break
        end
        % caso o criterio de convergencia nao tenha sido satisfeito, 
        % continuar com a proxima iteracao
        ypre=yh;
        
        % incrementando o contador de iteracoes
        it=it+1;
        % efetuando a reflexao
        xr=(1+ar)*xbar-ar*xh;
        % calculando a funcao no ponto refletido
        k=10.^xr;
        [p,qj]=calc_pres(t,n,phi,k',ct,h,muo,rw,q);
        yr=fobj(p,qj,pref,qref,flag);
        
        % entrando na arvore de decisao (fig. 6 de Luersen)
        if (yr<yl)
            % efetuando a expansao
            xe=ae*xr+(1-ae)*xbar;
            % desfazendo a transformacao logaritmica
            k=10.^xe;
            % calculando a funcao no ponto expandido
            [p,qj]=calc_pres(t,n,phi,k',ct,h,muo,rw,q);
            ye=fobj(p,qj,pref,qref,flag);
            % se ye < yh, entao a expansao foi bem sucedida
            if (ye < yl)
                % nesse caso, substituir xh por xe
                x(:,jh)=xe;
            % caso contrario, a expansao falhou
            else
                % nesse caso, substituir xh por xr
                x(:,jh)=xr;
            end
        % seguindo a arvore de decisao para o caso onde a reflexao falhou
        else
            % verificando se a reflexao gerou um novo maximo
            if (yr < yh)
            % se a reflexao NAO gerar novo maximo, realizar contracao com xr
                xc=ac*xr+(1-ac)*xbar;
            % se a reflexao gerar um novo maximo, realizar contracao com xh
            else
                xc=ac*xh+(1-ac)*xbar;
            end
            % desfazendo a transformacao logaritmica
            k=10.^xc;
            % calculando a funcao no ponto contraido
            [p,qj]=calc_pres(t,n,phi,k',ct,h,muo,rw,q);
            yc=fobj(p,qj,pref,qref,flag);
            % se yc < yh, entao a contracao foi bem sucedida
            if (yc < yh)
            % nesse caso, substituir xh por xc
                x(:,jh)=xc;
            else
            % caso a contracao tenha falhado, efetuar o encolhimento
                for jj=1:n+1
                    if (jj ~= jl)
                        x(:,jj)=as*(x(:,jj)+yl);
                    end
                end
            end
        end
    end
    % removendo as entradas vazias da matriz de xk e do vetor de residuos
    if (it < imax)
        xk(:,it+1:end)=[];
        res(it+1:end)=[];
    end
    res=real(res);
    % parando a contagem do tempo que leva para o algoritmo convergir
    t=toc;
end

%% funcao que calcula a pressao no poco e a vazao em cada camada
function [pwf, qj]=calc_pres(t, nlay, phi, k, ct, h, muo, rw, q)
    % obtendo o tamanho do vetor de tempo
    dim=length(t);
    % inicializando o vetor de pressoes com zeros
    pwf=zeros(dim,1);
    % inicializando a matriz de vazoes com zeros
    qj=zeros(dim,nlay);
    
    % constantes para ajuste de unidades
    alphat=0.0003484;
    alphap=19.03;
    
    % calculando a permeabilidade equivalente
    keq=sum(k.*h)/sum(h);
    
    for i=1:dim
        % calculando o valor de beta
        beta=4*alphat*t(i)/exp(.57722)/rw/rw/phi/ct;

        % calculando a pressao pela aproximacao logaritmica da linha fonte
        pwf(i)=alphap*q*muo/2/keq/sum(h)*(log(beta)+log(keq));
        
        % calculando as vazoes em cada camada
        for j=1:nlay
            omega=k(j)*h(j)/alphap/muo;
            qj(i,j)=pwf(i)*omega*2/(log(k(j))+log(beta));
        end
    end
end

%% funcao objetivo a ser minimizada
function[y]=fobj(p,qj,pref,qref,flag)
    % determinando o numero de camadas
    nlay=size(qref,2);
    % calculando o erro relativo dos dados de pressao
    y=norm(p-pref)/norm(pref);
    % se o flag for igual a 2, somar o erro relativo dos dados de vazao 
    if flag==2
        for j=1:nlay
            % para cada camada, somando o erro relativo dos dados de vazao
            y=y+norm(qj(:,j)-qref(:,j))/norm(qref(:,j));
        end
    end 
end

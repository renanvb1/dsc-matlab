% =========================================================================
%   Criado por Renan Vieira em 04/21.
%   Esse script aplica o algoritmo de Nelder-Mead para encontrar o
%   minimo local de uma funcao analitica definida
% =========================================================================

% limpando o console, a memoria e fechando as janelas graficas
clc; clear; close all

%% definindo alguns parametros de entrada
% definindo o flag da funcao analitica
flag=7;

% definindo a dimensao do problema (numero de parametros a serem estimados)
n=2;
% definindo o limite inferior em x1
x1=-2*ones(n,1);
% definindo o limite superior em x1
x2=-x1;

% definindo a tolerancia
tol=1e-10;
% definindo o numero maximo de iteracoes
nmax=150;

%% calculando o minimo usando o algoritmo de Nelder-Mead

[xest, xk, t, res]=NM_search(n,x1,x2,tol,nmax,flag); 

% exibindo os resultados
fprintf(' <strong> Algoritmo de Nelder-Mead </strong> \n')
% exibindo o tempo decorrido para executar o algoritmo
fprintf('Tempo decorrido: %.4f      No iteracoes: %d\n',t,length(res))
% exibindo as coordenadas do "chute inicial"
fprintf('Chute inicial:\n')
for ii=1:n
    fprintf('x%d: %.2f	',ii,xk(ii,1))
end
% exibindo o valor da funcao avaliada no "chute inicial"
fprintf('\n Valor da funcao em Xini: %.3f\n',f(xk(:,1),flag))
% exibindo as coordenadas do ponto estimado
fprintf('Ponto estimado:\n')
for ii=1:n
    fprintf('x%d: %.4f	',ii,xest(ii))
end
% exibindo o valor da funcao avaliada no ponto estimado
fprintf('\n Valor da funcao em Xest: %.6f\n',f(xest,flag))

%% exibindo alguns graficos
% calculando o valor da funcao a cada iteracao
zest=zeros(1,size(xk,2));
for jj=1:size(xk,2)
    zest(jj)=f(xk(:,jj),flag);
end

% plotando o residuo em cada iteracao
figure
yyaxis left
semilogy(res,'Linewidth',1.0)
yline(tol,'k--')
ylabel('Valor do residuo')
yyaxis right
semilogy(zest,'Linewidth',1.0)
ylabel('Funcao avaliada em xk')
xlabel('iterações')

%% para dimensao 2, plotar a superficie definida pela funcao analitica
if n==2
    
    % definindo o valor de dx e dy para plotar a superficie
    delta=0.1;
    % definindo a malha de referencia para plotar a superficie
    [mx,my]=meshgrid(x1:delta:x2,x1:delta:x2);
    
    % inicializando o valor de referencia z com zeros
    zref=zeros(size(mx));
    % para cada entrada da malha, calcular o valor de z
    for ii=1:size(mx,1)
        for jj=1:size(mx,2)
            zref(ii,jj)=f([mx(ii,jj) my(ii,jj)],flag);
        end
    end
    
    % plotando os dados
    figure
    s=surf(mx,my,zref);
    s.EdgeColor='none';
    colorbar
    hold on
    plot3(xk(1,:),xk(2,:),zest,'k*-','LineWidth',0.2,'MarkerSize',10)
    xlabel('Eixo x')
    ylabel('Eixo y')
end

%% limpando algumas variaveis auxiliares
clear ii jj delta zref mx my s

%% algoritmo de Nelder-Mead
function [xbar, xk, t, res]=NM_search(n,xlow,xhigh,tol,imax,flag)
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
    
    % entrando no loop do algoritmo
    while true
        % avaliando o valor da funcao em cada vertice do simplex inicial
        for jj=1:n+1
            y(jj)=f(x(:,jj),flag);
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
        % calculando o centroide de todos os pontos exceto xh
        for jj=1:n+1
            % somando todos os vertices
            xbar=xbar+x(:,jj);
        end
        % calculando a media sem xh
        xbar=(xbar-xh)/n;
        
        % atualizando o residuo
        ybar=(sum(y)-yh)/n;
        res(it)=sqrt(sum((y-ybar).^2)/n);
        
        % incluindo a estimativa atual na matriz de iteracoes
        xk(:,it)=xbar;
        
        % verificando se o numero maximo de iteracoes foi excedido
        % ou se o residuo esta menor que a tolerancia estabelecida
        if (it >= imax || res(it) < tol)
            break
        end
        % caso o criterio de convergencia nao tenha sido satisfeito, 
        % continuar com a proxima iteracao
        
        % incrementando o contador de iteracoes
        it=it+1;
        % efetuando a reflexao
        xr=(1+ar)*xbar-ar*xh;
        % calculando a funcao no ponto refletido
        yr=f(xr,flag);
        
        % entrando na arvore de decisao (fig. 6 de Luersen)
        if (yr<yl)
            % efetuando a expansao
            xe=ae*xr+(1-ae)*xbar;
            % calculando a funcao no ponto expandido
            ye=f(xe,flag);
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
            % calculando a funcao no ponto contraido
            yc=f(xc,flag);
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
    % parando a contagem do tempo que leva para o algoritmo convergir
    t=toc;
end

%% funcao que se deseja minimizar
function [y]=f(x,flag)
    switch flag
        case 1 % paraboloide centrado na origem
            y=sum(x.^2);
        case 2 % semiesfera de raio r e centro no ponto (0, 0, r)
            r=3;
            y=-sqrt(r*r-sum(x.^2))+r;
        case 3 % gaussiana de media 0 e desvio padrao sd
            sd=0.75;
            m=0;
            y=exp(-0.5*sum((x-m).^2)/sd/sd);
            y=-y/sqrt(2*pi*sd*sd)+0.6;
        case 4 % soma de 2 gaussianas de sd 0.5 e medias -0.75 e 0.75
            sd=.5;
            m1=0.75;
            m2=-m1;
            y=exp(-0.5*sum((x-m1).^2)/sd/sd);
            y=y+exp(-0.5*sum((x-m2).^2)/sd/sd);
            y=-y/sqrt(2*pi*sd*sd)+m1+0.1;
            y=y/10;
        case 5 % soma de cossenos
            arg=x/2*pi*1.5;
            y=sum(cos(arg));
        case 6 % soma do modulo de todas as coordenadas
            y=sum(abs(x));
        case 7 % produto do modulo de todas as coordenadas
            y=prod(abs(x));
        case 8 % produto do modulo de todas as coords + modulo de x1
            y=prod(abs(x))+sum(abs(x));
        otherwise % para flag nao especificado, retornar 0
            y=0;
    end
end


% =========================================================================
%   Created by Renan Vieira on 01/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   Esse script representa uma 2a tentativa de implementar a solucao
%   para testes de injetividade em reservatorios multicamadas com 
%   pocos horizontais multirramificados (condutividade infinita no poco)
% ====================================================== ===================

%% limpando o console e a memoria, e fechando as janelas graficas
clear; 
clc; close all
tic
%% sessao de definicao do vetor de tempo e das variaveis globais ==========
% definindo as constantes de conversao de unidades
global alphap alphat
% constante de conversao de unidades de pressao
alphap=19.03;
% constante de conversao de unidades de tempo
alphat=0.0003484;

% definido o flag que indica se há apenas injecao ou se há injecao e falloff
% Se o flag for igual a 1, há apenas injecao. Se o flag for igual a 2, há 
% injecao e falloff (na interface do STRITA, esse flag representa o numero
% de linhas na tabela de vazoes. Portanto, sera definido pelo usuario. OBS:
% por enquanto, a formulacao para poco horizontal permite que haja apenas
% 1 ou 2 periodos, ao contrario do caso com poco vertical, em que pode
% haver n periodos distintos. Alem disso, caso haja um segundo periodo, a
% vazao durante o segundo periodo sera obrigatoriamente zero, pois o modelo
% ainda nao consegue lidar com multiplas vazoes).
fper=2;

% definindo o flag que indica se o fluxo é monofasico ou bifasico (1 =
% monofasico; 2 = bifasico; No STRITA, o fluxo sempre sera bifasico,
% entao esse flag sera sempre igual a 2)
flap=2;

% definindo o tempo inicial
% (OBS: na interface atual do STIRTA, nao existe um campo para definicao do
% tempo inicial. Sugiro nao mexer nesse valor de t0 que eu defini abaixo)
t0=9.6e-8;
% definindo o tempo de fechamento do poco (definido pelo usuario)
% (No STRITA, o tempo de fechamento do poco vai ser o tempo de injecao 
% inserido na 1a linha da tabela de vazoes)
tp=96;
% definindo o numero de passos em cada periodo (Aqui eu defini
% arbitrariamente o numero de passos de tempo, e tambem determinei
% arbitrariamente que, nos casos em que ha periodo de injecao e de falloff,
% os 2 periodos tem o mesmo numero de passos de tempo. No STRITA, o numero 
% de passos em cada periodo é definido pelo usuario, na tabela de vazoes)
a=round(log10(tp/t0));
dim=10*a+1;

% calculando o vetor de tempo para o periodo de injecao
t=fill_time_new(t0,tp,dim);
% se houver periodo de falloff, criar os passos de tempo para o falloff
if fper==2
    % definido o tempo total do teste (aqui, eu defini arbitrariamente que
    % o tempo total de teste seria igual a 2 vezes o tempo de fechamento do
    % poco. Contudo, no STRITA, o tempo total de teste sera calculado com 
    % base na tabela de vazoes, a partir dos dados fornecidos pelo usuario
    tend=2*tp;
    % chamando a funcao para criar um vetor auxiliar com os passos de tempo do falloff
    taux=fill_time_new(t0,tend-tp,dim);
    taux=taux+tp;
    % acrescentando os passos de tempo do falloff no vetor de tempo
    t=[t; taux];
    % limpando a variavel com o vetor de tempo auxiliar
    clear taux
end

%% sessao de inicializacao dos parametros =================================
% essa sessao consiste simplesmente na inicializacao dos parametros de
% entrada (os dados que o usuario vai dar como input na interface)

% definindo os parametros que nao dependem da camada ======================
% definindo a vazao de injecao (definida na 1a linha da tabela de vazoes 
% definida pelo usuario. Por enquanto, o modelo considera apenas 2
% periodos: um onde a vazao é nao nula e outro onde a vazao é
% obrigatoriamente zero. Esse segundo periodo pode nao ser calculado, caso
% o usuario insira apenas uma linha na tabela de vazoes - no codigo, essa
% condicao é representada pela variavel fper)
qinj=1000;
% definindo o raio do poco (em m)
rw=0.0762;
% definindo o valor da viscosidade da agua 
miw=0.52;
% definindo o fator volume de formacao da agua (em m^3/ STD m^3)
bw=1;
% definindo o fator volume de formacao do oleo (em m^3/ STD m^3)
bo=1;
% definindo a pressao inicial (em kgf/cm^2)
p0=250;

%% inicializando os parametros de cada camada =============================
% definido o numero de camadas (será definido na interface do STRITA)
nlayers=2;
% (ao longo do codigo e em todas as funcoes, cada coluna representa uma camada)

% inicializando os parametros com um vetor de uns (no matlab, essa inicializacao 
% promove um ganho computacional. Os valores de fato sao definidos a partir
% da linha 134)
% permeabilidade na direcao x (no poco HORIZONTAL, as permeabilidades podem
% ser diferentes em cada direcao. Por isso, temos kx, ky e kz)
kx=ones(1,nlayers);
% permeabilidade na direcao y
ky=kx;          % parametro exclusivo para poco HORIZONTAL
% permeabilidade na direcao z
kz=kx;          % parametro exclusivo para poco HORIZONTAL
% altura da camada
h=kx;
% comprimento do poco na camada j
len=kx;        % parametro exclusivo para poco HORIZONTAL
% raio de skin
rskin=kx;
% permeabilidade de skin
kskin=kx;
% porosidade
phi=kx;
% definido a viscosidade do oleo (definido pelo usuario)
mio=kx;
% compressibilidade do oleo
co=kx;
% compressibilidade da rocha
cr=kx;
% compressibilidade da agua
cw=kx;
% o parametro dz abaixo representa a menor distancia entre o poco e alguma
% fronteira vertical (pode ser a fronteira superior ou a fronteira inferior
% dependendo da posicao do poco. Para o metodo, tanto faz qual delas)
dz=kx;        % parametro exclusivo para poco HORIZONTAL

%% sessao de definicao dos parametros de entrada===========================
% permeabilidades na direcao x (em mD) (definidos pelo usuario na interface)
kx(1)=300;     
kx(2)=220; 
% permeabilidades na direcao y (em mD) (definidos pelo usuario na interface)
% (a principio, a permeabilidade na direcao x pode ser diferente da
% permeabilidade na direcao y. Mas nesse exemplo, foi considerado que elas
% sao iguais)
ky=kx;
% permeabilidades na direcao z (em mD) (definidos pelo usuario na interface)
% (a principio, a permeabilidade na direcao z pode ser diferente da
% permeabilidade nas direcoes x e y. A razao entre kx e kz tambem pode ser
% diferente entre uma camada e outra. Nesse exemplo, eu defini
% arbitrariamente que os valores de kx, ky e kz fossem iguais.)
kz=kx/10;
% alturas das camadas (em m) (definidos pelo usuario na interface)
h(1)=20;       
h(2)=15;      
% comprimento do poco em cada camada (em m) (definidos pelo usuario na interface)
len(1)=350;    
len(2)=500;  
% porosidades (definidos pelo usuario na interface)
phi(1)=0.23;   
phi(2)=0.18;  
% raios de skin (em m) (definidos pelo usuario na interface)
rskin=0.56*rskin;
% permeabilidades de skin (em mD) (definidos pelo usuario na interface)
kskin(1)=13;
kskin(2)=2.5;  
% viscosidade do oleo (em cP) (definda pelo usuario na interface)
mio=3.1*mio;
% distancia entre o poco e uma fronteira vertical (definidos pelo usuario)
dz(1)=10;    dz(2)=7.5; 
% compressibilidade da rocha (definda pelo usuario na interface)
cr=1e-5*cr;
% compressibilidade do oleo (definda pelo usuario na interface)
co=1.0e-4*co;
% compressibilidade da agua (definda pelo usuario na interface)
cw=4.06e-5*cw;

%% sessao que busca os dados de permeabilidade relativa e mobilidade
% inicializando o flag que indica quais dados de perm relativa devem ser
% usados (o modelo permite que cada camada tenha uma curva de perm relativa
% diferente, mas nesse exemplo estou considerando que as curvas sao iguais
% em todas as camadas. Por isso o flag nao muda com a camada)
flag=0;
% no ambiente do STRITA, esses dados sao definidos pelo usuario, mas aqui
% eu busco os dados usando a funcao get_krel.
for jj=1:nlayers
    % obs: as matrizes com os dados de sw, krw e kro nao foram
    % inicializadas previamente, pois suas dimensoes dependem da quantidade
    % de pontos na curva de perm relativa. Por isso, nao foi feita a
    % prealocacao com zeros ou uns, como nos parametros anteriores
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(flag);
end

%% calculando as mobilidades extremais, a mobilidade total e o fluxo frac.
% A PARTIR DAQUI, TODOS OS DADOS DE ENTRADA DA INTERFACE DEVEM TER SIDO LIDOS

% DETALHE IMPORTANTE: embora as curvas de permeabilidade relativa possam
% ser diferentes em cada camada (como comentado acima), eu considerei que,
% em todas as camadas, o numero de entradas na curva de perm rel sera o
% mesmo. Teoricamente, cada camada pode ter uma quantidade de pontos na
% curva de perm rel diferente (e a interface atual do STRITA permite fazer
% isso). Mas por uma questao de praticidade, eu considerei que as curvas de
% perm rel vao ter o mesmo tamanho em todas as camadas. Do jeito que esta o
% codigo hoje, caso o usuario entre com curvas de perm rel de tamanhos
% diferentes em cada camada, o codigo vai dar erro de inconsistencia e nao
% vai conseguir fazer os calculos corretamente

% calculando as mobilidades extremais (esses valores sao definidos a partir 
% da tabela de permeabilidades relativas e das viscosidades.)
lohat=kro(1,:)./mio;
lwhat=krw(end,:)./miw;
% inicializando as matrizes de fluxo fracionario e mobilidade total (cada coluna representa uma camada)
lambdat=zeros(length(kro(:,1)),nlayers);
dfw=lambdat;
% calculando as matrizes de fluxo fracionario e mobilidade total (cada coluna representa uma camada)
for jj=1:nlayers
    [lambdat(:,jj),dfw(:,jj)]=fill_data(length(sw(:,jj)),krw(:,jj),kro(:,jj),sw(:,jj),mio(jj),miw);
end

%% sessao que calcula os parametros do reservatorio definidos a partir dos dados de entrada
% calculando a permeabilidade media de cada camada
% (esse parametro é definido como a raiz cubica de kx*ky*kz)
kj=(kx.*ky.*kz).^(1/3.0);
% calculando a compressibilidade total (ct = cr + cw*swi + co*(1-sor))
ct=cr+sw(1,:).*cw+sw(end,:).*co;
% calculando a difusividade hidraulica em cada camada (eta = k/(phi*mio*ct))
% (esse parametro é calculado a partir dos valores de kj, kro, ct, phi e mio)
etaj=kj.*kro(1)./phi./mio./ct;
% calculando os fatores de skin (S = (k/kskin - 1)*ln(rskin/rw))
Sj=(sqrt(kx.*kz)./kskin-1).*log(rskin./rw);

%% sessao onde sera chamada a funcao que calcula a pressao e as vazoes
% calculando os dados de pressao e vazao
[pwf, dpo, dpl, qj]=HW_press_new(t,tp,fper,flap,kj,kx,ky,kz,h,len,dz,...
    kskin,rskin,Sj,phi,qinj,etaj,dfw,ct,lambdat,lohat,lwhat,nlayers,rw,bo,bw);
% na chamada da funcao acima, os "..." apenas representam que ha mais
% argumentos de entrada para a mesma funcao na linha debaixo

% calculando a derivada da pressao
[dpwf]=compute_derivative(t,pwf,tp,fper);

%% plotando os resultados
% daqui pra baixo é apenas visualizacao de resultados

% definindo os limites dos eixos nos graficos
a=round(log10(t(1))); a=max([10^(a) 1e-3]); 
b=floor(log10(tp))+1; b=10^b;
c=31;
d=log10(mean(pwf(c:c+10)));
d=min([10^floor(d-1) .1]);
e=floor(log10(max(pwf)))+1; e=10^e;

% espessura da linha no grafico
linw=1.5;
% tamanho da fonte
fonts=14;
% plotando os dados durante o periodo de injecao
figure
loglog(t(c:dim),pwf(c:dim),'-k', t(c:dim),dpwf(c:dim),'-b','LineWidth',linw);
legend({'Pwf ', 'dPwf '}, 'Location','NorthWest', 'fontsize',fonts*.85)
axis([a,b,d,e])
title('Pressure and Pressure Derivative Profile During Injection')
xlabel('t (h)', 'fontsize',fonts)
ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
grid on
if fper==2
    % plotando os dados durante o periodo de injecao
    figure
    loglog(t(dim+c:end)-tp,pwf(dim)-pwf(dim+c:end),'-k', ...
        t(dim+c:end)-tp,dpwf(dim+c:end),'-b','LineWidth',linw);
    legend({'Pwf ', 'dPwf '}, 'Location','NorthWest', 'fontsize',fonts*.85)
    axis([a,b,d,e])
    title('Pressure and Pressure Derivative Profile During Falloff')
    xlabel('\Delta t (h)', 'fontsize',fonts)
    ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
    grid on
end

% plotando os dados de vazao
leg=strings(1,nlayers);
for jj=1:nlayers
    leg(jj)=strcat('Layer ',{' '},num2str(jj),{' '});
end
figure
semilogx(t(c:dim),qj(c:dim,:),'LineWidth',linw*1.1);
legend(leg, 'Location','West', 'fontsize',fonts*.85)
xlim([a,b])
title('Layer Flow-Rate Profiles')
xlabel('\Delta t (h)', 'fontsize',fonts)
ylabel('q (m³/d)', 'fontsize',fonts)
grid on

%% deletando algumas variaveis que nao sao mais necessarias
clear alphap alphat a b bo bw c co cr ct cw d dfw dim e etaj flag flap
clear fonts fper jj kj kro krw ky leg linw nlayers p0 sw t0 tend tp
toc
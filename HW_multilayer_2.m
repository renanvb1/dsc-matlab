% ==============================================================================================================
%   Created by Renan Vieira on 01/21.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   Esse script representa uma 2a tentativa de implementar a solucao
%   para testes de injetividade em reservatorios multicamadas com 
%   pocos horizontais multirramificados (condutividade infinita no poco)
%   ainda usando como base o modelo proposto no "analytical formulation"
% ==============================================================================================================

clear; 
clc; close all     % limpando o console e a memoria, e fechando as janelas graficas=============================

%% sessao de definicao do vetor de tempo e das variaveis globais ==========
% definindo as constantes de conversao de unidades
global alphap alphat
% constante de conversao de unidades de pressao
alphap=19.03;
% constante de conversao de unidades de tempo
alphat=0.0003484;

% definido o flag que indica se há apenas injecao ou se há injecao e falloff
% se o flag for igual a 1, há apenas injecao. Se o flag for igual a 2, há 
% injecao e falloff (definido a partir da tabela de vazoes definida pelo 
% usuario. OBS: por enquanto, a formulacao para poco horizontal permite que
% haja apenas 1 ou 2 periodos, ao contrario do caso com poco vertical, em
% que pode haver n periodos distintos. Alem disso, caso haja um segundo
% periodo, a vazao durante o segundo periodo sera obrigatoriamente zero,
% pois o modelo ainda nao consegue lidar com multiplas vazoes).
flagperiodo=2;

% flag para alternar o modo de inic. dos parametros (1 = txt; 2 = matlab)
flagi=1;
% definindo o flag que indica se o fluxo é monofasico ou bifasico (1 =
% monofasico; 2 = bifasico; No STRITA, o fluxo sempre sera bifasico,
% entao esse flag sera sempre igual a 2)
flap=2;
tp=96;
% definindo o tempo inicial
% (na interface atual do strita, nao existe um campo para o tempo inicial,
% entao eu defini um valor que pequeno o bastante para enxergar
% caracteristicas da pressao e da derivada no curto tempo)
t0=9.6e-8;
% definindo o tempo de fechamento do poco (definido pelo usuario)
% (No STRITA, o tempo de fechamento do poco vai ser o tempo de injecao 
% inserido na 1a linha da tabela de vazoes)
% tp=96;
a=round(log10(tp/t0));
% definindo o numero de passos em cada periodo (Aqui eu defini
% arbitrariamente que cada periodo tem o mesmo numero de passos de tempo. No 
% STRITA, o numero de passos sera definido pelo usuario, na tabela de vazoes)
dim=10*a+1;
% calculando o vetor de tempo para o periodo de injecao
t=fill_time_new(t0,tp,dim);
% se houver periodo de falloff, criar os passos de tempo para o falloff
if flagperiodo==2
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

% definindo o fator volume de formacao da agua (em m^3/ STD m^3)
bw=1;
% definindo o fator volume de formacao do oleo (em m^3/ STD m^3)
bo=1;
% definindo a pressao inicial (em kgf/cm^2)
p0=250;

% lendo os dados de entrada
% definindo o valor da viscosidade da agua (no STRITA, esse valor sera
% fornecido pelo usuario na interface)
miw=0.52;

if flagi==1
%% lendo os parametros de entrada no arquivo txt (no STRITA, todos os dados
% que aqui estao sendo lidos no arquivo txt devem ser fornecidos pelo
% usuario na interface)
% definindo o nome do arquivo
filename='props.txt';
% lendo o arquivo txt
fid=importdata(filename);
% usando uma funcao do matlab para salvar apenas os valores numericos no txt
fid=fid.data;

% determinando o numero de camadas a partir do arquivo txt (o numero de
% camadas nesse caso vai ser igual ao numero de linhas no arquivo txt)
nlayers=length(fid(:,1));
% lendo o raio do poco
rw=fid(1,10);
% lendo a vazao
qinj=fid(1,11);
Sj=[5 7 6 7];

% lendo as porosidades das camadas
phi=fid(:,1)';
% ao ler a porosidade, o simbolo ' serve apenas para transpor o vetor. Nao
% tem problema em armazenar os dados como vetores coluna, ao inves de
% vetores linha. Mas num ponto posterior do codigo, sera mais conveniente
% assumir que cada coluna representa uma camada. Entao, para ser
% consistente, eu optei por transpor esses vetores apenas para que aqui
% tambem cada coluna represente uma camada
% 
% lendo as permeabilidades na direcao x
kx=fid(:,2)';
% o modelo permite que as permeabilidades na direcao y sejam diferentes das
% permeabilidades na direcao x. Mas nesse caso, por conveniencia, assumi
% que kx = ky (no STRITA, a principio kx e ky sao fornecidos separadamente)
ky=kx;
% lendo as permeabilidades na direcao z
kz=fid(:,3)';
% lendo os raios da regiao de skin
rskin=fid(:,4)';
% lendo as permeabilidades na regiao de skin
kskin=fid(:,5)';
% lendo o comprimento do poco em cada camada
len=fid(:,6)';
% lendo as espessuras das camadas
h=fid(:,7)';
% definindo a menor distancia entre o poco e uma fronteira vertical (no
% STRITA, esse parametro sera fornecido pelo usuario)
dz=h./2; dz(2)=5.0;
% compressibilidade da rocha (definidos pelo usuario na interface)
% ct(1)=1.22e-4; ct(2)=ct(1);     ct(3)=ct(1);
ct=1.22e-4*ones(1,nlayers);

% lendo as viscosidades do oleo
mio=fid(:,8)';

% buscando os flags de permeabilidade relativa
flag=fid(:,9);
% (no STRITA, esses dados sao fornecidos pelo usuario na definicao da 
% camada. Aqui, eu defini um flag para cada camada que indica quais sao os 
% dados de perm rel corretos e usei a funcao abaixo para buscar os dados)


%% sessao de inicializacao dos parametros =================================
% essa sessao consiste simplesmente na inicializacao dos parametros de
% entrada (os dados que o usuario vai dar como input na interface)
else
% definindo os parametros que nao dependem da camada ======================
% definindo a vazao de injecao (definida na 1a linha da tabela de vazoes 
% definida pelo usuario. Por enquanto, o modelo considera apenas 2
% periodos: um onde a vazao é nao nula e outro onde a vazao é
% obrigatoriamente zero. Esse segundo periodo pode nao ser calculado, caso
% o usuario insira apenas uma linha na tabela de vazoes)
qinj=1500;
% definindo o raio do poco (em m)
rw=0.0762;

%% inicializando os parametros de cada camada =============================
% definido o numero de camadas (será definido na interface do STRITA)
nlayers=2;
% (ao longo do codigo e em todas as funcoes, cada coluna representa uma camada)
% inicializando os parametros com um vetor de uns (no matlab, essa inicializacao 
% promove um ganho computacional. Os valores de fato sao definidos a partir
% da linha 119)
% permeabilidade na direcao x (no poco HORIZONTAL, as permeabilidades podem
% ser diferentes em cada direcao)
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
% compressibilidade total
ct=kx;
% o parametro dz abaixo representa a menor distancia entre o poco e alguma
% fronteira vertical (pode ser a fronteira superior ou a fronteira inferior
% dependendo da posicao do poco. Para o metodo, tanto faz qual delas)
dz=kx;        % parametro exclusivo para poco HORIZONTAL

%% sessao de definicao dos parametros de entrada===========================
% permeabilidades na direcao x (em mD) (definidos pelo usuario na interface)
kx(1)=500;     
kx(2)=1450; 
% permeabilidades na direcao y (em mD) (definidos pelo usuario na interface)
% (a principio, a permeabilidade na direcao x pode ser diferente da
% permeabilidade na direcao y. Mas nesse exemplo, foi considerado que elas
% sao iguais)
ky=kx;
% permeabilidades na direcao z (em mD) (definidos pelo usuario na interface)
kz=(kx/2);
% alturas das camadas (em m) (definidos pelo usuario na interface)
h(1)=20;       
h(2)=15;      
% comprimento do poco em cada camada (em m) (definidos pelo usuario na interface)
len(1)=100;    
len(2)=120;  
% porosidades (definidos pelo usuario na interface)
phi(1)=0.23;   
phi(2)=0.27;     
% phi(3)=0.14;

% raios de skin (em m) (definidos pelo usuario na interface)
rskin=0.35*rskin;
% permeabilidades de skin (em mD) (definidos pelo usuario na interface)
kskin(1)=200;
kskin(2)=200;  
% viscosidade do oleo (em cP) (definda pelo usuario na interface)
mio=5.7*mio;
% compressibilidade da rocha (definidos pelo usuario na interface)
% ct(1)=1.22e-4; ct(2)=ct(1);     ct(3)=ct(1);
ct=1.22e-4*ct;
% distancia entre o poco e uma fronteira vertical (definidos pelo usuario)
dz(1)=10.0;    dz(2)=5.0; 
end
%% sessao que busca os dados de permeabilidade relativa e mobilidade
% inicializando o flag que indica quais dados de perm relativa devem ser
% usados (o modelo permite que cada camada tenha uma curva de perm relativa
% diferente, mas nesse exemplo estou considerando que as curvas sao iguais
% em todas as camadas, por isso o flag nao muda com a camada)
if flagi==2
flag=1;
for jj=1:nlayers
    % obs: as matrizes com os dados de sw, krw e kro nao foram
    % inicializadas previamente, pois suas dimensoes dependem da quantidade
    % de pontos na curva de perm relativa. Por isso, nao foi feita a
    % prealocacao com zeros ou uns, como nos parametros anteriores
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(flag);
end
else
% no ambiente do STRITA, esses dados sao definidos pelo usuario, mas aqui
% eu busco os dados usando a funcao get_krel.
for jj=1:nlayers
    % obs: as matrizes com os dados de sw, krw e kro nao foram
    % inicializadas previamente, pois suas dimensoes dependem da quantidade
    % de pontos na curva de perm relativa. Por isso, nao foi feita a
    % prealocacao com zeros ou uns, como nos parametros anteriores
    [sw(:,jj),krw(:,jj),kro(:,jj)]=get_krel(flag(jj));
end
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
% calculando a difusividade hidraulica em cada camada
% (esse parametro é calculado a partir dos valores de kj, kro, ct, phi e mio)
etaj=kj.*kro(1)./phi./mio./ct;
% if flagi==2
% calculando os fatores de skin
Sj=(sqrt(kx.*kz)./kskin-1).*log(rskin./rw);
% end

%% sessao onde sera chamada a funcao que calcula a pressao e as vazoes
% calculando os dados de pressao e vazao
[pwf, dpo, dpl, qj]=HW_press_new(t,tp,flagperiodo,flap,kj,...
    kx,ky,kz,h,len,dz,kskin,rskin,Sj,phi,qinj,etaj,dfw,ct,...
    lambdat,lohat,lwhat,nlayers,rw,bo,bw);
% na chamada da funcao acima, os "..." apenas representam que ha mais
% argumentos de entrada para a mesma funcao na linha debaixo

% calculando a derivada da pressao
[dpwf]=compute_derivative(t,pwf,tp,flagperiodo);

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
seq=sum(kzx.*len.*(Sj+1*sani))/kzxeq/sum(len);

% expoente ej necessario para o calculo do skin total (lrf)
ej=kxy.*h./kzx./len;
% Sgj necessario para o calculo do skin total (lrf)
sgj=len.*((2*pi*rw./h).^ej);     sgj=log(4*rw./sgj);
% skin total em cada camada (que ocorre durante o lrf)
stj=sgj+kxy.*h./kzx./len.*Sj;
% skin total equivalente
steq=sum(kxy.*h.*stj)/kxyeq/sum(h);

% patamares teoricos de derivada
merf=alphap*qinj/2/lohatm/sum(len)/kzxeq;
mlrf=alphap*qinj/2/lohatm/sum(h)/kxyeq;

%% calculando as aproximacoes logaritmicas de curto e longo tempo
c=10*(a-6)+1;
% definindo o vetor de tempo para o curto tempo
terf=t(c:c+21);
% definindo o vetor de tempo para o longo tempo
tlrf=t(dim-11:dim);

% definido o argumento do ln no curto tempo
argerf=4*kzxeq*lohatm*alphat/exp(.57722)/phicteq/rw/rw;
% definido o argumento do ln no longo tempo
arglrf=4*kxyeq*lohatm*alphat/exp(.57722)/phicteq/rw/rw;

% calculando as aproximacoes logaritmicas de curto e longo tempo
perf=merf*(log(argerf.*terf)+2*seq);
plrf=mlrf*(log(arglrf.*tlrf)+2*steq);

%% plotando os resultados
% daqui pra baixo é apenas visualizacao de resultados

% definindo os limites dos eixos nos graficos
a=round(log10(t(1))); a=max([10^(a) 1e-3]); 
b=floor(log10(tp))+1; b=10^b;

d=min(min(abs(pwf(c:dim))),min(abs(dpwf(c:dim)))); 
d=10^floor(log10(d));
e=floor(log10(max(pwf)))+1; e=10^e;

linw=1.5;
fonts=14;
% plotando os dados durante o periodo de injecao
figure
loglog(t(c:dim),pwf(c:dim),'-k', t(c:dim),dpwf(c:dim),'-b','LineWidth',linw);
hold on
plot(terf,perf,'gv', tlrf,plrf,'rv','LineWidth',linw)
hold off
yline(merf,'k--')
yline(mlrf,'k:','LineWidth',1.2*linw)
legend({'Pwf ', 'dPwf '}, 'Location','NorthWest', 'fontsize',fonts*.85)
axis([a,b,d,e])
title('Pressure and Pressure Derivative Profile During Injection')
xlabel('t (h)', 'fontsize',fonts)
ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
grid on
if flagperiodo==2
    % plotando os dados durante o periodo de injecao
    figure
    loglog(t(dim+c:end)-tp,pwf(dim)-pwf(dim+c:end),'-k', ...
        t(dim+c:end)-tp,dpwf(dim+c:end),'-b','LineWidth',linw);
%     hold on
%     plot(terf,perf,'gv', tlrf,plrf,'rv','LineWidth',linw)
%     hold off
    yline(merf,'k--')
    yline(mlrf,'k:','LineWidth',1.2*linw)
    legend({'Pwf ', 'dPwf '}, 'Location','NorthWest', 'fontsize',fonts*.85)
    axis([a,b,d,e])
    title('Pressure and Pressure Derivative Profile During Falloff')
    xlabel('\Delta t (h)', 'fontsize',fonts)
    ylabel('\Delta P (kgf/cm²)', 'fontsize',fonts)
    grid on
end

p1=pwf(c:dim);
p2=pwf(dim+c:end);
qjinj=qj(c:dim,:);
tinj=t(c:dim);

%% deletando algumas variaveis que nao sao mais necessarias
clear alphap alphat bo bw ctj dim etaj flagperiodo kj p0
clear ct dfw kro krw ky nlayers rfzx xf rfxy qj t 
clear a b c d e jj flag fonts linw argerf arglrf fid filename
clear kzx kxy kzxeq kxyeq ej sgj stj perf plrf phicteq terf tlrf
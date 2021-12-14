% funcao que retorna os dados de permeabilidade relativa
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


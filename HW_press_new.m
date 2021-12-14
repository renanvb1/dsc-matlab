% funcao que calcula os dados de pressao e vazoes
function[pwf, dpo, dpl, qj]=HW_press_new(t,tp,flaper,flagp,kj,kx,ky,kz,h,len,dz,...
    kskin,rskin,Sj,phi,qinj,eta,dfw,ct,lt,lohat,lwhat,nlayers,rw,bo,bw)
    % buscando a constante global de conversao de unidades alphap
    global alphap
        
    %% inicializando os parametros de saida
    % inicializando os dados de pressao com zeros
    pwf=zeros(size(t));
    % o pwf representa a pressao medida no poco. Ele é definido como a soma
    % de 2 termos: o deltaPo (representado pelo vetor dpo) e o deltaPlambda
    % (representado pelo vetor dpl). A cada passo de tempo, esses termos
    % sao calculados separadamente e depois somados
    dpo=pwf;
    dpo2=dpo;
    dpl=pwf;
    calcdpo=1;
    % na definicao da matriz de vazao, cada coluna representa uma camada
    qj=zeros(length(t),nlayers);
    
    %% calculos preliminares
    
    % determinando qual o fator volume de formacao deve ser empregado nas contas
    if flaper==1
        bf=bo;
    else
        bf=bw;
    end
    
    % calculando a permeabilidade equivalente do reservatorio no plano xy
    keqxy=sqrt(kx.*ky);
    keqxy=sum(keqxy.*h)/sum(h);
    % calculando o produto porosidade-compressibilidade equivalente
    phict=sum(phi.*ct.*h)./sum(h);
    
    % o modelo para poco horizontal apresenta uma hipotese diferente do que
    % ocorre na realidade. Para compensar essa diferenca no modelo, a
    % pressao verdadeira deve ser calculada a partir da pressao media que
    % sera obtida pelo modelo. Entao, precisamos definir o numero de
    % particoes ao longo do poco que serao feitas para calcularmos a
    % pressao media. Definido o numero de particoes no poco: (pode diminuir
    % esse numero para testar a implementacao, mas depois deve voltar para
    % o valor original) 15 monof
    np=15;
    
    % por hipotese, o poco esta posicionado paralelo ao eixo y. A pressao
    % sera calculada num ponto de coordenadas deltax, deltay, deltaz. Para
    % simplificar a notacao, o ponto sera representado por dx, dy, dz. Na
    % direcao y, precisamos definir um conjunto np de pontos para poder
    % calcular a pressao media ao longo do poco (que esta paralelo a y).
    % Tambem precisamos definir um conjunto de pontos para cada camada pois
    % existe uma ramificacao do poco em cada camada
    
    % definindo as posicoes dx e dz (cada coluna representa uma camada)
    dx=rw*ones(1,nlayers);
    deltaz=dx;
    % inicializando as posicoes dy
    dy=zeros(np,nlayers);
    % definindo as posicoes dy (cada coluna representa uma camada)
    if np~=1
        for jj=1:nlayers
            % calculando o incremento que sera dado a cada dy
            dyaux=len(jj)/(np-1);
            for ii=1:np
                dy(ii,jj)=dyaux*(ii-1);
            end
        end
    else
        for jj=1:nlayers
            dy(1,jj)=len(jj)/2;
        end
    end
    
    % inicializando a matriz com os coeficientes Aj (cada coluna representa uma camada)
    Aj=qj;
    
    % calculando o comprimento do vetor de tempo
    dim=length(t);
    % encontrando o ponto de tempo em que t = tp
    if flaper==2
        % se houver periodo de falloff, buscando o tempo em que t = tp
        flagt=find(t>tp);
        flagt=flagt(1);
    else
        % se houver apenas o periodo de injecao, definindo flagt como a
        % ultima posicao do vetor de tempo
        flagt=length(t);
    end
    
    % o modelo proposto calcula a pressao pelo metodo de fontes e
    % sumidouros. Por esse metodo, precisamos primeiro calcular a variacao
    % de pressao instantanea provocada pela fonte. Essa variacao
    % instantanea depende de t, do dy e da camada (pois cada camada tem uma
    % ramificacao do poco). Por isso, as variacoes instantaneas das fontes
    % foram armazenadas numa "matriz 3D", que pode ser construida no matlab
    % a coordenada i representa o tempo, a coordenada j representa a camada
    % e a coordenada k representa o dy
    
    % inicializando a matriz 3D com as variacoes instantaneas
    inst=zeros(dim, nlayers, np);
    % calculando as entradas da matriz 3D com as variacoes instantaneas
    for n=1:np
        for ii=1:length(t)
            inst(ii,:,n)=HW_inst_sourc(t(ii),eta,kx,ky,kz,kj,...
                len,dx,dy(n,:),deltaz,h,dz);
        end
    end
    % Para calcular a pressao, deve ser resolvido 1 sistema linear A*x = b,
    % onde o vetor b é constante mas os coeficientes da matriz A mudam a
    % cada passo de tempo. Por isso, o sistema precisa ser recalculado a
    % cada passo de tempo. A matriz A é definida através de relacoes entre
    % os coeficientes sigma. Esses coeficientes sao especificos de cada
    % camada (portanto, representado por um vetor) e mudam com o tempo. Com
    % isso, a matriz A tambem muda no tempo.
    % inicializando os vetores sigma, b e a matriz A
    % o vetor sigma foi inicializado com zeros e sera calculado a cada
    % passo de tempo por meio de uma integral. Ele tem dimensao 1 X nlayers
    sigma=zeros(1,nlayers);
    % os coeficientes da matriz A serao atualizados a cada passo de tempo
    % todos os elementos na primeira linha da matriz A sao constantes e
    % iguais a 1. Os demais elementos dependem das relacoes entre os
    % coeficientes sigma, que serao calculados a cada tempo. Por isso, na
    % inicializacao da matriz A, a primeira linha foi preenchida com uns e
    % as demais linhas foram preenchidas com zeros. A dimensao da matriz é
    % nlayers X nlayers
    A=[ones(1,nlayers);
        zeros(nlayers-1,nlayers)];
    % inicializando o vetor b com zeros. Ele é constante, e apenas o
    % primeiro termo é nao-nulo (ver linha abaixo). O vetor b tem dimensao
    % nlayers X 1
    b=zeros(nlayers,1);
    % o primeiro elemento do vetor b deve ser igual a qinj
    b(1)=qinj;
    
    %% calculo da pressao e das vazoes
    fc=1;
    if fc==2
        % definindo a posicao auxiliar ao longo do poco (necessaria para as contas)
        naux=round(.75*np);
    end
    % entrando no loop de tempo durante o periodo de injecao
    for ii=1:flagt
        % calculando o coeficiente sigma no passo de tempo atual
        if fc==2
            for j=1:nlayers
                sigma(j)=HW_comp_sigma(ii,inst(:,j,naux),t,len(j),phi(j),ct(j),bf);
            end
        else
            for j=1:nlayers
                sigma(j)=0;
                % para cada camada, calculando o coeficiente sigma
                for n=1:np
                    sigma(j)=sigma(j)+HW_comp_sigma(ii,inst(:,j,n),t,len(j),phi(j),ct(j),bf);
                end
            end
            sigma=sigma/np;
        end
        % atualizando os elementos da matriz A (lembrando que os elementos
        % da primeira linha de A devem ser sempre constantes e iguais a 1)
        for j=1:nlayers-1
            A(j+1,j)=sigma(j)+alphap*bf*Sj(j)/sqrt(kx(j)*kz(j))/lohat(j)/len(j);
%             A(j+1,j)=sigma(j);
            if j+1<=nlayers
                A(j+1,j+1)=-sigma(j+1)-alphap*bf*Sj(j+1)/sqrt(kx(j+1)*kz(j+1))...
                    /lohat(j+1)/len(j+1);
%                 A(j+1,j+1)=-sigma(j+1);
            end
        end
        % resolvendo o sistema linear no passo de tempo atual para
        % determinar as vazoes nesse passo de tempo
        qjaux=A\b;
        % transpondo o vetor de vazoes e salvando na 1a linha da matriz
        qj(ii,:)=qjaux';
%         % calculando as demais entradas do vetor b
%         if nlayers>1
%             for jaux=2:nlayers
%                 % definindo um indice de tempo auxiliar para evitar acessar
%                 % indices negativos no vetor de tempo
%                 ia=max(ii-1,1);
%                 % calculando as entradas do vetor
%                 b(jaux)=qj(ia,jaux)*sigma(jaux)-qj(ia,jaux-1)*sigma(jaux-1);
%             end
%         end
%         % resolvendo o sistema linear para determinar as vazoes no tempo ii
%         qjaux=A\b;
%         % transpondo o vetor de vazoes e salvando os dados na matriz de vazoes
%         qj(ii,:)=qjaux';
        
        %% calculando o termo dpl e os coeficientes Aj no tempo ii
        if flagp==2
            [dpl(ii), Aj(ii,:)]=HW_Aj_Rj(ii, rw, qinj, t, qj, h, kx, ky, kz,...
                keqxy, phi, phict, len, kskin, rskin, lohat, lwhat, lt, dfw, tp, dz, 1, bf);
        end
        %% calculando o termo deltapskin
        dpskin=alphap*qj(ii,1)*bf*Sj(1)/sqrt(kx(1)*kz(1))/len(1)/lohat(1);
        if nlayers>1 && calcdpo==1
            dpski2=alphap*qj(ii,2)*bf*Sj(2)/sqrt(kx(2)*kz(2))/len(2)/lohat(2);
        end
        
        %% calculando o termo deltapo
        for n=1:np
        % somando o valor de deltapo calculado em cada subparticao do poco
            dpo(ii)=dpo(ii)+HW_deltapo_new(ii,inst(:,1,n),t,len(1),phi(1),...
                ct(1),qj(:,1),bf);
            if nlayers>1 && calcdpo==1
                dpo2(ii)=dpo2(ii)+HW_deltapo_new(ii,inst(:,2,n),t,len(2),...
                    phi(2),ct(2),qj(:,2),bf);
            end
        end
        % calculando deltapo a partir da media ao longo do poco
        dpo(ii)=dpo(ii)/np+dpskin;
        if nlayers>1 && calcdpo==1
            dpo2(ii)=dpo2(ii)/np+dpski2;
        end
        
        %% calculando a pressao no poco como a soma dos 2 termos
        pwf(ii)=dpo(ii)+dpl(ii);
    end
    if nlayers>1 && calcdpo==1
    dpski2=(dpo(31:flagt)-dpo2(31:flagt))./dpo(31:flagt)*100;
    disp(dpski2)
    figure
    yyaxis left
    loglog(t(31:flagt),dpo(31:flagt),'k',t(31:flagt),dpo2(31:flagt),'b')
    hold on
    xlim([1e-3 1e2])
    ylim([1e-1 1e1])
    yyaxis right
    plot(t(31:flagt),abs(dpski2))
    ylim([0 20])
    end
    
    %% se houver falloff, definindo as matrizes auxiliares
    if flagt~=length(t)
        % definindo a matriz auxiliar com as vazoes em t=tp
        qjfall=qj(flagt,:).*ones(dim-flagt,nlayers);
%         sigaux=zeros(1,nlayers);
%         for jj=1:nlayers
%             qjaux(:,jj)=qj(flagt,jj)*qjaux(:,jj);
%         end
        % definindo o vetor de tempo auxiliar
        taux=t(flagt+1:end)-tp;
        % inicializando a matriz com as fontes instantaneas auxiliares
        instaux=zeros(dim-flagt, nlayers, np);
        % calculando as fontes instantaneas auxiliares
        for n=1:np
            for ii=1:length(taux)
                instaux(ii,:,n)=HW_inst_sourc(taux(ii),eta,kx,ky,kz,kj,...
                    len,dx,dy(n,:),deltaz,h,dz);
            end
        end
    end
    
    %% entrando no loop de tempo durante o falloff
    for ii=flagt+1:dim
        % inicializando as variaveis auxiliares para determinar dpo
        dpa=0.0;
        dpb=0.0;
        % modificando provisoriamente a vazao para efetuar os calculos
        qj(ii,:)=qj(flagt,:);
        
%         for j=1:nlayers
%             sigma(j)=0;
%             % para cada camada, calculando o coeficiente sigma
%             for n=1:np
%                 sigma(j)=sigma(j)+HW_comp_sigma(ii,inst(:,j,n),t,len(j),phi(j),ct(j),bf);
%             end
%         end
%         sigma=sigma/np;
%         % atualizando os elementos da matriz A (lembrando que os elementos
%         % da primeira linha de A devem ser sempre constantes e iguais a 1)
%         for j=1:nlayers-1
%             A(j+1,j)=sigma(j)+alphap*bf*Sj(j)/sqrt(kx(j)*kz(j))/lohat(j)/len(j);
% %             A(j+1,j)=sigma(j);
%             if j+1<=nlayers
%                 A(j+1,j+1)=-sigma(j+1)-alphap*bf*Sj(j+1)/sqrt(kx(j+1)*kz(j+1))...
%                     /lohat(j+1)/len(j+1);
% %                 A(j+1,j+1)=-sigma(j+1);
%             end
%         end
        % resolvendo o sistema linear no passo de tempo atual para
%         % determinar as vazoes nesse passo de tempo
%         qjaux=A\b;
%         % transpondo o vetor de vazoes e salvando na 1a linha da matriz
%         qj(ii,:)=qjaux';
%         disp(ii)
        % somando o valor de deltapo em cada subparticao do poco
        for n=1:np
            % definindo a contribuicao a partir do tempo inicial
            dpa=dpa+HW_deltapo_new(ii,inst(:,1,n),t,len(1),phi(1),ct(1),...
                qj(:,1),bf);
            dpb=dpb+HW_deltapo_new(ii-flagt,instaux(:,1,n),taux,len(1),phi(1),...
                ct(1),qjfall(:,1),bf);
        end
        % calculando deltapo usando superposicao e a media ao longo do poco
        dpo(ii)=dpa/np-dpb/np;
        % calculando deltaplambda durante o falloff
        if flagp==2
            [dpl(ii), Aj(ii,:)]=HW_Aj_Rj(ii, rw, qinj, t, qj, h, kx, ky, kz,...
                keqxy, phi, phict, len, kskin, rskin, lohat, lwhat, lt, dfw, tp, dz, 2, bf);
        end
        
        %% calculando a pressao no poco como a soma dos 2 termos
        pwf(ii)=dpo(ii)+dpl(ii);
    end
    % corrigindo a matriz de vazoes
    if flagt~=length(t)
        qj(flagt+1:end,:)=zeros(length(t)-flagt,nlayers);%qj(flagt+1:end,:)-qjaux;
    end
end


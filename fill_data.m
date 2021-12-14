% funcao que recebe os dados de permeabilidade relativa e viscosidade para
% calcular as curvas de mobilidade total e derivada do fluxo fracionario
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


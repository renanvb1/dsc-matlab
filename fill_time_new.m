% funcao que cria o vetor de tempo
function [t]=fill_time_new(t0,tend,npassos)
    % os inputs da funcao sao o tempo inicial, o tempo final e o no de passos
    
    % inicializando o vetor de tempo com zeros
    t=zeros(npassos,1);
    % a funcao preenche o vetor de tempo conforme uma progressao geometrica
    % calculando a razao da PG de acordo com os tempos inicial, final e o no de passos
    tim=(tend/t0)^(1/(npassos-1));
    
    % definindo o tempo inicial na primeira posicao do vetor
    t(1)=t0;
    % preenchendo as outras posicoes do vetor conforme a PG
    for ii=2:npassos
        t(ii)=t(ii-1)*tim;
    end
end
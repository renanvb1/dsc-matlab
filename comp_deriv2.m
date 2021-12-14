function [dpwf]=comp_deriv2(t,pwf,tp,flagperiodo)
% Funcao que calcula a derivada da pressao com relacao ao ln do tempo
    
    % definindo o tamanho do vetor
    dim=length(pwf);
    % initializando o vetor de derivada com zeros
    dpwf=zeros(size(pwf));
    
    % definindo o espacamento
    gap=4;
    
    % calculando a derivada nos 1os pontos
    dpwf(1)=(pwf(2)-pwf(1))/log(t(2)/t(1));
    for ii=2:gap
        % calculando os 2 termos necessarios para obter a derivada de Bourdet
        a = (pwf(ii + 1) - pwf(ii)) / log(t(ii + 1) / t(ii))*log(t(ii) / t(ii - 1));
        b = (pwf(ii) - pwf(ii - 1)) / log(t(ii) / t(ii - 1))*log(t(ii + 1) / t(ii));
        % Derivada de Bourdet disponivel no pdf 74 do arquivo "notas cap 3"
        dpwf(ii)=(a+b)/log(t(ii + 1) / t(ii - 1));
    end
%     dpwf(1)=(pwf(2)-pwf(1))/log(t(2)/t(1));
    
    % encontrando o ponto de tempo em que t = tp
    if flagperiodo==2
        % se houver periodo de falloff, buscando o tempo em que t = tp
%         tflag=find(abs(t-tp)<1e-10);
        tflag=find(t>tp);
        tflag=tflag(1);
        % para evitar problemas de aproximacao numerica, a busca tolera um erro de 1e-12
        
        % para todos os pontos com pelo menos 2 vizinhos, calculando a
        % derivada segundo Bourdet
        for ii=gap+1:tflag-gap-1
            % calculando os 2 termos necessarios para obter a derivada de Bourdet
            a = (pwf(ii + gap) - pwf(ii)) / log(t(ii + gap) / t(ii))*log(t(ii) / t(ii - gap));
            b = (pwf(ii) - pwf(ii - gap)) / log(t(ii) / t(ii - gap))*log(t(ii + gap) / t(ii));
            % Derivada de Bourdet disponivel no pdf 74 do arquivo "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + gap) / t(ii - gap));
        end
        % calculando a derivada nos ultimos pontos
        for ii=tflag-gap:tflag
            dpwf(ii)=(pwf(ii)-pwf(ii-1))/log(t(ii)/t(ii-1));
        end
        
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
        for ii=gap+1:tflag-gap-1
            % calculando os 2 termos necessarios para obter a derivada de Bourdet
            a = (pwf(ii + gap) - pwf(ii)) / log(t(ii + gap) / t(ii))*log(t(ii) / t(ii - gap));
            b = (pwf(ii) - pwf(ii - gap)) / log(t(ii) / t(ii - gap))*log(t(ii + gap) / t(ii));
            % Derivada de Bourdet disponivel no pdf 74 do arquivo "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + gap) / t(ii - gap));
        end
        % calculando a derivada nos ultimos pontos
        for ii=tflag-gap:tflag-1
            % calculando os 2 termos necessarios para obter a derivada de Bourdet
            a = (pwf(ii + 1) - pwf(ii)) / log(t(ii + 1) / t(ii))*log(t(ii) / t(ii - 1));
            b = (pwf(ii) - pwf(ii - 1)) / log(t(ii) / t(ii - 1))*log(t(ii + 1) / t(ii));
            % Derivada de Bourdet disponivel no pdf 74 do arquivo "notas cap 3"
            dpwf(ii)=(a+b)/log(t(ii + 1) / t(ii - 1));
        end
        dpwf(tflag)=(pwf(tflag)-pwf(tflag-1))/log(t(tflag)/t(tflag-1));
    end
    
end
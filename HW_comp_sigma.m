% function that computes deltapo at a given time using the sources and sinks method
function [sigma] = HW_comp_sigma(index, sources, t, len, phi, ct, b)
    flag=1;
    global alphap
    % calculando a integral no 1o passo de tempo
    sigma=sources(1)*t(1);
    % incrementando a integral
    for ii=2:index
        % calculando numericamente a integral usando a regra do trapezio
        sigma=sigma+(sources(ii-1,:)+sources(ii,:)).*(t(ii)-t(ii-1))/2;
    end
    % multiplicando a integral por uma constante
    if flag==2
        sigma=sigma*alphap*b/(24*len*phi*ct);
    else
        sigma=sigma*b/(24*len*phi*ct);
    end
end
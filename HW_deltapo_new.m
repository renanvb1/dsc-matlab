% funcao que calcula o deltapo em um dado tempo 
function [dpo] = HW_deltapo_new(index, inst, t, len, phi, ct, qj, b)
    flag=1;
    global alphap
    % calculando a integral no 1o passo de tempo
    dpo=inst(1,:).*t(1).*qj(1,:);
    % incrementando a integral nos demais tempos
    for ii=2:index
        % calculando a integral usando a regra do trapezio
        dpo=dpo+(inst(ii-1,:)*qj(ii-1,:)+inst(ii,:)*qj(ii,:))*(t(ii)-t(ii-1))/2;
%         dpo=dpo+(inst(ii,:)*qj(ii,:))*(t(ii)-t(ii-1));
    end
    % multiplicando o deltapo por uma constante
    if flag==2
        dpo=dpo*alphap*b/(24*len.*phi.*ct);
    else
        dpo=dpo*b/(24*len.*phi.*ct);
    end
end
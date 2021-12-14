function [dpl, aj] = HW_Aj_Rj(ii, rw, qinj, t, qj, h, kx, ky, kz, keq,...
    phi, phict, len, kskin, rskin, lohat, lwhat, lt, dfw, tp, dz, flagper, bf)
% funcao que calcula o termo deltaplambda para pocos horizontais em um
% determinado tempo. Essa funcao NAO foi modificada em relacao a versao
% anterior, houve alteracao apenas nos comentarios

    %% inicializacoes
    % buscando a constante global de conversao de unidades alphap
    global alphap
    % para cada camada, sera necessario calcular 3 parcelas: uma no plano
    % zx, uma na direcao x e uma no plano xy
    
    % inicializando a parcela no plano zx como zero
    dpzx=zeros(size(len));
    % inicializando a parcela na regiao de skin como zero
    dpzx_s=dpzx;
    % inicializando a parcela na direcao x como zero
    dpx=dpzx;
    % inicializando a parcela no plano xy como zero
    dpxy=dpzx;
    
    % para o calculo de cada parcela, sera necessario calcular a frente de
    % avanco correspondente a cada direcao
    
    % inicializando a frente de avanco no plano zx com zeros
    rfzx=zeros(size(dfw));
    % inicializando a frente de avanco na direcao x com zeros
    xf=rfzx;
    % inicializando a frente de avanco no plano xy com zeros
    rfxy=rfzx;

    %% calculando o termo deltaplambda durante o periodo de injecao========
    if (flagper==1)
        %% calculando a frente de avanco em cada camada, em cada direcao
        for jj=1:length(len)
            % calculando a frente de avanco no plano zx em cada camada
            rfzx(:,jj)=compute_rf(ii, rw, t, qj(:,jj), len(jj), phi(jj), dfw(:,jj));
            % calculando a frente de avanco na direcao x em cada camada
            xf(:,jj)=compute_xf(ii, rw, t, qj(:,jj), h(jj), len(jj), phi(jj), dfw(:,jj));
            % calculando a frente de avanco no plano xy em cada camada
            rfxy(:,jj)=compute_rf(ii, 0.0, t, qj(:,jj), h(jj), phi(jj), dfw(:,jj));
        end
        %% calculando deltapzx segundo eq 48 do arquivo "09-SPE-84957"
        
        % definindo o limite superior da 1a integral na eq 48 do arquivo "09-SPE-84957"
        a=min(dz,rfzx(1,:));
        % para cada camada, calculando deltapzx
        for j=1:length(len)
            % definindo o primeiro valor do indice
            ind=length(rfzx(:,j));
            % calculando a 1a integral na eq 48 do arquivo "09-SPE-84957"
            while (rfzx(ind,j)<a(j) && ind>0)
                raux=min(rfzx(ind-1,j), a(j));
                % incrementando a integral usando a regra do trapezio
                dpzx(j)=dpzx(j)+log(raux/rfzx(ind,j))*(lohat(j)/lt(ind,j)...
                    -1+lohat(j)/lt(ind-1,j)-1)/2;
                % atualizando o indice
                ind=ind-1;
            end
            
            % calculando a 2a integral na eq 48 do arquivo "09-SPE-84957"
            ind=length(rfzx(:,j));
            while (rfzx(ind,j)<rskin(j)&&ind>1)
                % checando se o proximo raio esta alem da regiao de skin
                if (rfzx(ind-1,j)>rskin(j))
                    % interpolando a mobilidade total
                    lamskin=lt(ind,j)+(lt(ind-1,j) - lt(ind,j))*...
                        (rskin(j) - rfzx(ind,j)) / (rfzx(ind-1,j) - rfzx(ind,j));
                    % incrementando a integral auxiliar usando a mobilidade interpolada
                    dpzx_s(j)=dpzx_s(j)+log(rskin(j)/rfzx(ind,j))*...
                        (lohat(j)/lt(ind,j)-1+lohat(j)/lamskin-1)/2;
                else
                    % incrementando a integral auxiliar usando a regra do trapezio
                    dpzx_s(j)=dpzx_s(j)+log(rfzx(ind-1,j)/rfzx(ind,j))*(lohat(j)/lt(ind,j)-1+lohat(j)/lt(ind-1,j)-1)/2;
                end
                % atualizando o indice
                ind=ind-1;
            end
        end
        % multiplicando a 2a integral por (k/kskin-1)
        dpzx_s=dpzx_s.*(sqrt(kx.*kz)./kskin-1);
        % somando as duas integrais e multiplicando pela constante
        dpzx=(dpzx+dpzx_s).*h./len./sqrt(kx.*kz);
        
        %% calculando deltapx segundo o arquivo "09-SPE-84957"
        
        % calculando x1 segundo o padrao geometrico proposto por Deppe
        x1=pi*h/8;
        % calculando x3 segundo o padrao geometrico proposto por Deppe
        x3=pi*len/8;
        % definindo o limite superior da integral na eq 49 do arquivo "09-SPE-84957"
        b=min(x3, max([x1;xf]));
        % calculando x2 segundo a eq 61 do arquivo "10-SPE-90907"
        x2=pi*pi/8./(h./2./dz-1);
        x2=x2.*log(h./2./dz)+log(h./2./pi./dz./sin(pi*dz./h));
        x2=x2.*h/pi;
        x2=x2./(log(h./2./dz).*h./2./dz./(h./2./dz-1)-1);
        % para cada camada, calculando dpx
        for j=1:length(len)
            % buscando o indice que corresponde ao limite inferior da integ
            ind=find(xf(:,j)>=x1(j));
            ind=max(ind);
            % calculando deltapx segundo a eq 49 do arquivo "09-SPE-84957"
            if isempty(ind)==0
                while (xf(ind,j)<b(j) && ind>0)
                    % calculando espessura efetivamente inundada em x1 e x2
                    if (xf(ind,j)<x2(j) && ~isnan(x2(j)) )
                        h1=h(j)-(h(j)-2*dz(j))./(x2(j)-x1(j))*(x2(j)-xf(ind));
                    else
                        h1=h(j);
                    end
                    if (xf(ind-1,j)<x2(j) && ~isnan(x2(j)) )
                        h2=h(j)-(h(j)-2*dz(j))./(x2(j)-x1(j))*(x2(j)-xf(ind-1));
                    else
                        h2=h(j);
                    end
                    % incrementando a integral usando a regra do trapezio
                    dpx(j)=dpx(j)+((lohat(j)/lt(ind,j)-1)/h1+(lohat(j)/lt(ind-1,j)-1)/h2)/2*(xf(ind-1,j)-xf(ind,j));
                    % atualizando o indice
                    ind=ind-1;
                end
            end
        end
        % multiplicando a integral by pi/L
        dpx=dpx.*pi./len./kx;
        
        %% calculando o termo deltapxy segundo o arquivo "09-SPE-84957"
        
        % calculando o limite superior da integral
        c=max([len./2;rfxy]);
        % para cada camada, calculando dpxy
        for j=1:length(len)
            % encontrando o indice que corresponde ao limite inferior 
            ind=find(rfxy(:,j)>=len(j)/2);
            ind=max(ind);
            % calculando deltapx segundo a eq 49 do arquivo "09-SPE-84957"
            if isempty(ind)==0
                while (rfxy(ind,j)<c(j) && ind>0)
                    % incrementando a integral usando a regra do trapezio
                    dpxy(j)=dpxy(j)+(lohat(j)/lt(ind,j)-1+lohat(j)/lt(ind-1,j)-1)/2*log(rfzx(ind-1,j)./rfzx(ind,j));
                    % atualizando o indice
                    ind=ind-1;
                end
            end
        end
        dpxy=dpxy./h./sqrt(kx.*ky);
        
        %% calculando o coeficiente Aj em cada camada
        
        % definindo Aj como a soma das contribuicoes nas 3 direcoes(eq 47 do arquivo "09-SPE-84957")
        aj=dpzx+dpx+dpxy;
        % multiplicando Aj por alphap/(kj*lohat*hj)
        aj=aj*alphap*bf./lohat./h;
        
        % calculando deltaplambda
        dpl=qinj/sum(1./aj);
        
    %% calculando o deltaplambda durante o periodo de falloff
    else
        for jj=1:length(len)
            % calculando a frente de avanco no plano zx
            rfzx(:,jj)=compute_rf(ii, rw, t, qj(:,jj), len(jj), phi(jj), dfw(:,jj));
            % calculando a frente de avanco na direcao x
            xf(:,jj)=compute_xf(ii, 0.0, t, qj(:,jj), h(jj), len(jj), phi(jj), dfw(:,jj));
            % calculando a frente de avanco no plano xy
            rfxy(:,jj)=compute_rf(ii, 0.0, t, qj(:,jj), h(jj), phi(jj), dfw(:,jj));
        end
        % definindo o flag para a aproximacao da vazao total
        approx=1;
        
        %% calculando a 1a integral em Rj
        % definindo o limite superior da 1a integral em Rj
        a=min(dz,rfzx(1,:));
        rdim=rw;
        % para cada camada, calculando a 1a integral em Rj
        for j=1:length(len)
            % calculando a integral entre os 2 primeiros raios 
            % numero de subparticoes entre os 2 primeiros raios (pode diminuir esse
            % numero para testar a implementacao, mas depois deve voltar para o
            % valor original) 10
            n=5;
            raux=min(rfzx(end-1,j), a(j));
            % incremento de cada subparticao entre os 2 primeiros raios
            dr=(raux-rw)/n;
            % calculando a integral entre os 2 primeiros raios
            for r1=rw:dr:raux-dr
                r2=r1+dr;
                lt1=lt(end,j)+(r1-rw)*(lt(end-1,j)-lt(end,j))/(rfzx(end-1,j)-rw);
                lt2=lt(end,j)+(r2-rw)*(lt(end-1,j)-lt(end,j))/(rfzx(end-1,j)-rw);
                qos1=approxlog(ii, r1-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
                qos2=approxlog(ii, r2-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
                qts1=approxlog(ii, r1-rdim, t, tp, rw, qinj, keq, lt1, phict);
                qts2=approxlog(ii, r2-rdim, t, tp, rw, qinj, keq, lt2, phict);
                dpzx(j)=dpzx(j)+log(r2/r1)*(qts1/lt1-qos1/lohat(j)+qts2/lt2-qos2/lohat(j))/2;
            end
            % definindo o primeiro valor do indice
            ind=length(rfzx(:,j))-1;
            while(rfzx(ind,j)<a(j) && ind>0)
                % definindo o numero de particoes entre os raios
                n=2;
                raux=min(rfzx(ind-1,j), a(j));
                % calculando o incremento entre cada subparticao dos raios
                dr=(raux-rfzx(ind,j))/n;
                % calculando a integral
                for r1=rfzx(ind,j):dr:raux-dr
                    r2=r1+dr;
                    lt1=lt(ind,j)+(r1-rw)*(lt(ind-1,j)-lt(ind,j))/(rfzx(ind-1,j)-rw);
                    lt2=lt(ind,j)+(r2-rw)*(lt(ind-1,j)-lt(ind,j))/(rfzx(ind-1,j)-rw);
                    qos1=approxlog(ii, r1-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
                    qos2=approxlog(ii, r2-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
                    qts1=approxlog(ii, r1-rdim, t, tp, rw, qinj, keq, lt1, phict);
                    qts2=approxlog(ii, r2-rdim, t, tp, rw, qinj, keq, lt2, phict);
                    dpzx(j)=dpzx(j)+log(r2/r1)*(qts1/lt1-qos1/lohat(j)+qts2/lt2-qos2/lohat(j))/2;
                end
                % atualizando o indice
                ind=ind-1;
            end
            % checando se o raio de skin esta alem da zona de vazao transiente
            qosskin=approxlog(ii, rskin(j)-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
            % calculando a contribuicao da zona de skin para o Rj
            if (qosskin>1e-10)
                pos=length(rfzx(:,j))-length(find(rfzx(:,j)<rskin(j)))+1;
                % determinnando a mobilidade total mobility no raio de skin
                ltskin=lt(pos,j)+(rskin(j)-rfzx(pos,j))*...
                    (lt(pos-1,j)-lt(pos,j))/(rfzx(pos-1,j)-rfzx(pos,j));
                % definindo o numero de subparticoes entre o poco e o raio de skin
                n=2;
                % calculando a integral do poco ate o raio de skin
                dr=(rskin(j)-rw)/n;
                for r1=rw:dr:rskin(j)-dr
                    r2=r1+dr;
                    lt1=lt(end,j)+(r1-rfzx(pos,j))*(ltskin-lt(end,j))/(rfzx(pos-1,j)-rfzx(pos,j));
                    lt2=lt(end,j)+(r2-rfzx(pos,j))*(ltskin-lt(end,j))/(rfzx(pos-1,j)-rfzx(pos,j));
                    qos1=approxlog(ii, r1-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
                    qos2=approxlog(ii, r2-rdim, t, tp, rw, qinj, keq, lohat(j), phict);
                    qts1=approxlog(ii, r1-rdim, t, tp, rw, qinj, keq, lt1, phict);
                    qts2=approxlog(ii, r2-rdim, t, tp, rw, qinj, keq, lt2, phict);
                    dpzx_s(j)=dpzx_s(j)+log(r2/r1)*(qts1/lt1-qos1/lohat(j)+qts2/lt2-qos2/lohat(j))/2;
                end
            end
        end
        % multiplicando dpzx_s por (k/kskin-1)
        dpzx_s=dpzx_s.*(sqrt(kx.*kz)./kskin-1);
        % somando as duas integrais e multiplicando por 1/sqrt(kz*kx)
        dpzx=(dpzx+dpzx_s)./sqrt(kx.*kz);
        
        %% calculando a 2a integral em Rj
        % calculando x1 segundo o padrao geometrico proposto por Deppe
        x1=pi*dz/8;
        % calculando x3 segundo o padrao geometrico proposto por Deppe
        x3=pi*len/8;
        % calculando x2 segundo a eq 61 do arquivo "10-SPE-90907"
        x2=pi*pi/8./(h./2./dz-1);
        x2=x2.*log(h./2./dz)+log(h./2./pi./dz./sin(pi*dz./h));
        x2=x2.*h/pi;
        x2=x2./(log(h./2./dz).*h./2./dz./(h./2./dz-1)-1);
        % defindo o limite superior da integral na eq 49 do arquivo "09-SPE-84957"
        b=min(x3, max([x1;xf]));
        % para cada camada, calculando dpx
        for j=1:length(len)
            % buscando o indice que corresponde ao limite superior da integ
            ind=find(xf(:,j)>=x1(j));
            ind=max(ind);
            % calculando deltapx segundo a eq 49 do arquivo "09-SPE-84957"
            if isempty(ind)==0
                while (xf(ind,j)<b(j) && ind>0)
                    % calculando a espessura efetiva em x1 e x2
                    if (xf(ind,j)<x2(j) && ~isnan(x2(j)) )
                        h1=h(j)-(h(j)-2*dz(j))./(x2(j)-x1(j))*(x2(j)-xf(ind));
                    else
                        h1=h(j);
                    end
                    if (xf(ind-1,j)<x2(j) && isnan(x2(j))==0)
                        h2=h(j)-(h(j)-2*dz(j))./(x2(j)-x1(j))*(x2(j)-xf(ind-1));
                    else
                        h2=h(j);
                    end
                    % atualizando as vazoes total e de oleo
                    qos1=qf_linear(ii, xf(ind,j), t, tp, qinj, keq, lohat(j), phict);
                    qos2=qf_linear(ii, xf(ind-1,j), t, tp, qinj, keq, lohat(j), phict);
                    qts1=qos1;
                    qts2=qos2;
                    % incrementando a integral usando a regra do trapezio
                    dpx(j)=dpx(j)+((qts1/lt(ind,j)-qos1/lohat(j))/h1+(qts2/lt(ind-1,j)-qos2/lohat(j))/h2)/2*(xf(ind-1,j)-xf(ind,j));
                    % atualizando o indice
                    ind=ind-1;
                end
            end
        end
        % multiplicando a integral por pi/kx
        dpx=dpx.*pi./kx;
        
        %% calculando a 3a integral em Rj
        % definindo o limite superior da integral
        c=max([len./2;rfxy]);
        % para cada camada, calculando dpxy
        for j=1:length(len)
            % buscando o indice que corresponde ao limite superior da integ
            ind=find(rfxy(:,j)>=len(j)/2);
            ind=max(ind);
            % calculando deltapx segundo a eq 49 do arquivo "09-SPE-84957"
            if isempty(ind)==0
                while (rfxy(ind,j)<c(j) && ind>0)
                    % atualizando as vazoes total e de oleo
                qos1=approxlog(ii, rfzx(ind-1,j), t, tp, rw, qinj, keq, lohat(j), phict);
                qos2=approxlog(ii, rfzx(ind,j), t, tp, rw, qinj, keq, lohat(j), phict);
                qts1=calcqts(ii, rfzx(ind-1,j), approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
                qts2=calcqts(ii, rfzx(ind,j), approx, t, tp, rw, qinj, keq, lohat(j), phict, lwhat);
                    % incrementando a integral usando a regra do trapezio
                    dpxy(j)=dpxy(j)+(qts1/lt(ind,j)-qos1/lohat(j)+qts2/lt(ind-1,j)-qos2/lohat(j))/2*log(rfzx(ind-1,j)./rfzx(ind,j));
                    % atualizando o indice
                    ind=ind-1;
                end
            end
        end
        dpxy=dpxy.*len./h./sqrt(ky.*kx);
        
        %% calculando Rj como a soma dos 3 termos
        aj=dpzx+dpx+dpxy;
        % multiplicando Rj por alphap/L
        aj=aj.*alphap*bf./len;
        
        % calculando deltaPlambda durante o periodo de falloff
        dpl=1/sum(1./aj);
    end
end


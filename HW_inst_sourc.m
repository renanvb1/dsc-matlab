% funcao que calcula as variacoes instantaneas de pressao devido a
% existencia do poco. O modelo usa o metodo de fontes e sumidouros
function [sources] = HW_inst_sourc(t,eta,kx,ky,kz,k,lenj,dx,dy,deltaz,h,dz)
    % recuperando o valor da constante global de conversao de unidades
    global alphat
    
    % no calculo, sera necessario fazer um somatorio infinito
    % definindo o limite computacional para o infinito (pode diminuir esse
    % numero para testar a implementacao, mas depois deve voltar para o
    % valor original) 250
    infinity=50*5;
    
    % calculando as difusividades hidraulicas em cada direcao
    etax=eta.*kx./k;
    etay=eta.*ky./k;
    etaz=eta.*kz./k;
    
    % pelo modelo, a variacao instantanea de pressao pode ser calculada
    % pelo produto entre tres termos: deltapx, deltapy e deltapz, que serao
    % representados por dpx, dpy e dpz
    
    % calculando o vetor com os valores de dpx em cada tempo
    dpx=exp(-dx.*dx./4./etax./alphat./t);
    dpx=dpx/2.0./sqrt(pi.*etax.*alphat.*t);
    
    % calculando o vetor com os valores de dpy em cada tempo
    dpy=erf((lenj-dy)./2.0./sqrt(etay.*alphat.*t))-erf(-dy./2.0./sqrt(etay.*alphat.*t));
    dpy=dpy/2.0;
    
% % %     % inicializando dpz como zero
% % %     dpz=0.0;
% % %     % fazendo a soma dos infinitos termos usando a serie de cosenos
% % %     for n=1:infinity
% % %         dpz=dpz+exp(-n*n*pi*pi*etaz.*alphat.*t./h./h).*cos(n.*pi./2).*cos(n.*pi./2);
% % %     end
% % %     dpz=(dpz*2+1)./h;
    
% %     % inicializando dpz como zero
% %     dpz=0.0;
% %     % fazendo a soma dos infinitos termos
% %     for n=-infinity:infinity
% %         % atualizando o coeficiente do somatorio
% %         an=n*h;
% %         % incrementando o somatorio
% %         dpz=dpz+exp(-(deltaz-an).*(deltaz-an)./4.0./etaz./alphat./t);
% %     end
% %     dpz=dpz/2./sqrt(pi.*etaz.*alphat.*t);

%     % para o calculo de dpz, sera necessario fazer um somatorio infinito
%     % calculando o parametro auxiliar Dz
%     Dz=h-dz;
%     % inicializando o coeficiente do somatorio infinito como sendo zero
%     an=0.*Dz+0.*dz;
%     % acrescentando o 1o termo do somatorio
%     dpz=exp(-(deltaz-an).*(deltaz-an)./4./etaz./alphat./t);
%     % fazendo a soma dos termos com coeficiente negativo
%     for n=-1:-1:-infinity
%         % atualizando o coeficiente do somatorio
%         an=an-2*mod(n,2).*Dz-2*mod(n+1,2).*dz;
%         % incrementando o somatorio
%         dpz=dpz+exp(-(deltaz-an).*(deltaz-an)./4./etaz./alphat./t);
%     end
%     % zerando o coeficiente do somatorio
%     an=0.*Dz+0.*dz;
%     % fazendo a soma dos termos com coeficiente positivo
%     for n=1:infinity
%         % atualizando o coeficiente do somatorio
%         an=an+2*mod(n,2).*Dz+2*mod(n+1,2).*dz;
%         % incrementando o somatorio
%         dpz=dpz+exp(-(deltaz-an).*(deltaz-an)./4./etaz./alphat./t);
%     end
%     dpz=dpz/2./sqrt(pi.*etaz.*alphat.*t);
    
    % para o calculo de dpz, sera necessario fazer um somatorio infinito
    % calculando o parametro auxiliar Dz
    Dz=h-dz;
    % inicializando dpz com zero
    dpz=zeros(size(dpx));
    for n=-infinity:infinity
        % atualizando o parametro auxiliar do somatorio
        an=2*(n-1)*dz+2*n*Dz;
        % incrementando o somatorio
        dpz=dpz+exp(-(deltaz-an).*(deltaz-an)./4./etaz./alphat./t);
        % atualizando o parametro auxiliar do somatorio
        an=2*n*dz+2*n*Dz;
        % incrementando o somatorio
        dpz=dpz+exp(-(deltaz-an).*(deltaz-an)./4./etaz./alphat./t);
    end
    dpz=dpz/2./sqrt(pi.*etaz.*alphat.*t);
    
    % calculando a variacao de pressao instantanea como sendo o produto de
    % dpx, dpy e dpz
    sources=dpx.*dpy.*dpz;
end


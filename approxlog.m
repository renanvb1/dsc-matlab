% function that computes the flow-rate at a given point of the reservoir using the line-source logarithmic approximation
function [qd]=approxlog(index, ray, t, tp, rw, qinj, keq, lambdahat, phict)
    global alphat
    % % initializing the dimensionless flow-rate as zero
    % qd=0.0;
%     % calculating the dimensionless radius
%     rd=ray/rw;
%     % computing the dimensionless times
%     td=alphat*keq*lambdahat*t(index)/(phict*rw*rw);
%     deltatd=alphat*keq*lambdahat*(t(index)-tp)/(phict*rw*rw);
%     % calculating the dimensionless flow-rate with respect to td using the logarithmic approximation
%     qd=exp(-rd*rd /(4*td));
%     % during falloff, adjust the dimensionless flow-rate to account for deltatd's contribuition
%     if t(index)>tp
%         qd=qd-exp(-rd*rd/4/deltatd);
%     end
    % computing the first argument
    arg1=phict*ray*ray/4/alphat/keq/lambdahat/t(index);
    % computing the second argument
    arg2=arg1*t(index)/(t(index)-tp);
    % computing the dimensionless flow-rate
    qd=exp(-arg1);
    % during falloff, computing the flow-rate using superposition 
    if t(index)>tp
        qd=qd-exp(-arg2);
    end
    % dimensionalizing qd
    qd=qd*qinj;
end
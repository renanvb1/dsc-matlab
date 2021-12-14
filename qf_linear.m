% function that computes the flow-rate of a given phase for a linear flow
function [qf] = qf_linear(index, x, t, tp, qinj, k, lambdahat, phict)
    global alphat
    % computing the complementary error fuction argument
    arg=x*sqrt(phict)/2;
    arg=arg/sqrt(alphat*k*lambdahat*t(index));
    % computing the flow-rate according to eq 18 from file "10-SPE-90907"
    qf=erfc(arg);
    % during falloff, applying the superposition principle
    if (t(index)>tp)
        % updating the complementary error function argument
        arg=x*sqrt(phict)/2;
        arg=arg/sqrt(alphat*k*lambdahat*(t(index)-tp));
        % applying the superposition principle to compute the flow-rate
        qf=qf-erfc(arg);
    end
    % dimensionalizing the flow-rate
    qf=qinj*qf;
end
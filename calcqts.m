% function that computes the total flow-rate at a given radius during falloff period
function [qts]=calcqts(index, ray, approx, t, tp, rw, qinj, keq, lohatm, phict, lwhatm)
    % initializing qts as zero
    qts=0.0;
    if (approx==1)
        % according to the 1st aproximation, qts = qos
        qts=approxlog(index, ray, t, tp, rw, qinj, keq, lohatm, phict);
    end
    if (approx==2)
        % according to the 2nd approximation, qts = qws
        qts=approxlog(index, ray, t, tp, rw, qinj, keq, lwhatm, phict);
    end
    if (approx==3)
        % according to the 3rd approximation, qts is an averaged mean using water and oil properties @ their endpoint saturations
        qos=approxlog(index, ray, t, tp, rw, qinj, keq, lohatm, phict);
        qws=approxlog(index, ray, t, tp, rw, qinj, keq, lwhatm, phict);
        qts=(qos*lohatm+qws*lwhatm)/(lohatm+lwhatm);
    end    
end
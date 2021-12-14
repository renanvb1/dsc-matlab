function [rfxy] = compute_rf(ii, rw, t, q, h, phi, dfw)
%     function computes the waterfront radius in the x-y plane
    % computing the integral at the first timestep
    integral=q(1,:).*t(1);
    % numerically integrating (using the trapeze rule) the expression required to determine waterfront radius
    for u=1:ii-1
        integral=integral+(t(u+1)-t(u)).*(q(u+1,:)+q(u,:))./2;
    end
    % calculating the waterfront radii using eq. 54 from file "09-SPE-84957"
    rfxy=sqrt(integral.*dfw./24./pi./phi./h+rw*rw);
end
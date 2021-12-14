function [xf] = compute_xf(ii, rw, t, q, h, len, phi, dfw)
%     function computes the waterfront position in the x direction
    % computing the integral at the first timestep
    integral=q(1,:).*t(1);
    % numerically integrating (using the trapeze rule) the expression required to determine waterfront radius
    for u=1:ii-1
        integral=integral+(t(u+1)-t(u)).*(q(u+1,:)+q(u,:))./2;
    end
    % calculating the waterfront radii using eq. 54 from file "09-SPE-84957"
    xf=integral.*dfw./24./2./phi./h./len+rw;
end
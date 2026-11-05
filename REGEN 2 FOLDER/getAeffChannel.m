function A = getAeffChannel(x,geom,geom_channels)

    r   = geom.r(x);
    drdx  = geom.drdx(x);
    kap = geom.kappa(x);
    dxds = geom.dxds(x);

    % aprime = d(dx/ds)/dx
    aprime = geom.aprime(x);

    tw = geom_channels.tw(x);
    hc = geom_channels.hc(x);
    wc = geom_channels.wc(x);
    dtwdx = geom_channels.dtwdx(x);
    dhcdx = geom_channels.dhcdx(x);
    dwcdx = geom_channels.dwcdx(x);

    % dtw/ds
    dtwds = dtwdx.*dxds;

    % Channel angle
    rho_midchan = r + dxds.*(tw + hc./2);
    theta_c = wc./rho_midchan;

    drhodx_midchan = drdx + aprime .* (tw + hc./2) + dxds.* (dtwdx + dhcdx./2);
    
    dthetacdx = (dwcdx.*rho_midchan - wc.*drhodx_midchan) ./ rho_midchan.^2;

    dthetacds =  dthetacdx.*dxds;
 
   

    rho = @(zeta) r + dxds.*zeta.*tw;
    Jg  = @(zeta) 1 + kap.*zeta.*tw;

    integrand_c = @(eta,zeta) ...
        tw .* rho(zeta) .* theta_c .* Jg(zeta) ./ ...
        (Jg(zeta).^2 ...
        + (zeta.*dtwds).^2 ...
        + (rho(zeta).*eta.*dthetacds).^2);


    Aeff_c = integral2(integrand_c, -0.5,0.5, 0,1);

    A = Aeff_c;
end
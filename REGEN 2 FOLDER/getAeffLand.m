function A = getAeffLand(x,geom,geom_channels,theta_p)

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
    theta_l = theta_p - theta_c;

    drhodx_midchan = drdx + aprime .* (tw + hc./2) + dxds.* (dtwdx + dhcdx./2);
    
    dthetacdx = (dwcdx.*rho_midchan - wc.*drhodx_midchan) ./ rho_midchan.^2;
    dthetaldx =  -dthetacdx;
    dthetalds =  dthetaldx.*dxds;
   

    rho = @(zeta) r + dxds.*zeta.*tw;
    Jg  = @(zeta) 1 + kap.*zeta.*tw;

    integrand_l = @(eta,zeta) ...
        tw .* rho(zeta) .* theta_l .* Jg(zeta) ./ ...
        (Jg(zeta).^2 ...
        + (zeta.*dtwds).^2 ...
        + (rho(zeta).*eta.*dthetalds).^2);

    Aeff_l = integral2(integrand_l, -0.5,0.5, 0,1);


    A = Aeff_l;

end
function val = axialRgIntegrandChannel(x,eta,zeta,geom,gch)

    r     = geom.r(x);
    drdx  = geom.drdx(x);
    kap   = geom.kappa(x);
    a     = geom.dxds(x);

    tw    = gch.tw(x);
    hc    = gch.hc(x);
    wc    = gch.wc(x);

    dtwdx = gch.dtwdx(x);
    dhcdx = gch.dhcdx(x);
    dwcdx = gch.dwcdx(x);

    aprime = geom.aprime(x);

    dtwds = dtwdx .* a;

    rho_mid = r + a.*(tw + hc/2);

    theta = wc ./ rho_mid;

    drhomiddx = ...
        drdx ...
        + aprime.*(tw + hc/2) ...
        + a.*(dtwdx + dhcdx/2);

    dthetadx = ...
        (dwcdx.*rho_mid - wc.*drhomiddx) ...
        ./ rho_mid.^2;

    dthetads = dthetadx .* a;

    rho = r + a.*zeta.*tw;
    J   = 1 + kap.*zeta.*tw;

    f = ...
        tw .* rho .* theta .* J ./ ...
        (J.^2 ...
        + (zeta.*dtwds).^2 ...
        + (rho.*eta.*dthetads).^2);

    val = sqrt(1 + drdx.^2) ./ f;

end
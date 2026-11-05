function Tz = wallTempProfile(zeta, Thot, Tmid, Tcold)

    Tz = zeros(size(zeta));

    lower = zeta <= 0.5;
    upper = ~lower;

    Tz(lower) = ...
        Thot + 2*zeta(lower).*(Tmid-Thot);

    Tz(upper) = ...
        Tmid + 2*(zeta(upper)-0.5).*(Tcold-Tmid);

end
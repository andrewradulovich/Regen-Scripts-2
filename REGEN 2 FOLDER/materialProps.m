%% Material Properties

% inputs
% material options: 'al6061-ram2'

function mat = materialProps(material)

   
    if material == "al6061-ram2" 
        properties = readmatrix(pwd + "/Material Data/al6061_RAM2.xlsx");
        E_data = [properties(17:22, 8) properties(17:22,9)];
        yts_data = properties(2:12,1:2);
        elong_data = [properties(2:12,1) properties(2:12,4)];


        % CTE Curve x: Kelvin y: ppm *10^-6 / K
        piece_CTE.curvefit = @(x) (1.87098e-15).*x.^4 + (-1.78017e-12).*x.^3 + (4.3767e-10).*x.^2 + (2.69185e-8).*x + (7.10145e-6);
        piece_CTE.extrapol = 21.87e-6;
        piece_CTE.extrastart = 373.15;
        CTE = @(x) CTE_func(x,piece_CTE);

        % k curve x: Kelvin y: W/m K
        % I have data to 452K, I think its reasonable to let it extrapolate
        % to 550K. This is where AL6061-RAM2 thermal conductivity peaks,
        % and begins decreasing as the material reaches its melting
        % temperature. I will assume that at 550K there is an inflection
        % point and the thermal conductivity starts decreasing at the same
        % slope for the purposes of extrapolation. 

        piece_k.ascending = @(T) 0.16463.*T + 145.354;
        piece_k.inflection = 550;
        k_inflection = piece_k.ascending(piece_k.inflection);
        piece_k.decending = @(T) k_inflection - 0.16463.*(T-550);
        k = @(T) k_func(T,piece_k);

        % Conductivity antiderivative
        piece_K_ad.ascending = @(T) (0.16463/2).*(T.^2) + 145.354.*T;
        piece_K_ad.inflection = 550;
        K_inflection = piece_K_ad.ascending(piece_K_ad.inflection);
        piece_K_ad.decending = @(T) K_inflection + k_inflection.*(T-550) - (0.16463/2).*((T-550).^2);
        K_ad = @(T) K_ad_func(T,piece_K_ad);

        % Youngs Modulus [Pa]
        piece_E.interp = @(T) interp1(E_data(:,1), E_data(:,2), T, 'linear');
        piece_E.lowertempbound = E_data(1,1);
        piece_E.uppertempbound = E_data(end,1);
        piece_E.lowerbound = E_data(1,2);
        piece_E.upperbound = E_data(end,2);
        E = @(x) E_func(x,piece_E);

        % Elongation [fractional] 
        piece_elong.interp = @(T) interp1(elong_data(:,1), elong_data(:,2), T, 'linear');
        piece_elong.lowertempbound = elong_data(1,1);
        piece_elong.uppertempbound = elong_data(end,1);
        piece_elong.lowerbound = elong_data(1,2);
        piece_elong.upperbound = elong_data(end,2);
        elong = @(x) elong_func(x,piece_elong);

        % Yield Strength [Mpa]
        piece_yts.interp = @(T) interp1(yts_data(:,1), yts_data(:,2), T, 'linear');
        piece_yts.lowertempbound = yts_data(1,1);
        piece_yts.uppertempbound = yts_data(end,1);
        piece_yts.lowerbound = yts_data(1,2);
        piece_yts.upperbound = yts_data(end,2);
        yts = @(x) yts_func(x,piece_yts);

    elseif material == "inco-718"
        % properties = readmatrix(pwd + "/Material Data/Inconel718.xlsx");
        % k = properties(13:end,1:2); % thermal conductivity [W/m-K]
        % yts = properties(1:7,1:2);
        % elong = [properties(1:7,1) properties(1:7,5)];
        % E = [properties(1:6, 9) properties(1:6,10)]; %youngs modulus
        % CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm]

        %% new **REDO***
        properties = readmatrix(pwd + "/Material Data/Inconel718.xlsx");
        E_data = [properties(1:6, 9) properties(1:6,10)];
        yts_data = properties(1:7,1:2);
        elong_data = [properties(1:7,1) properties(1:7,5)];


        % CTE Curve x: Kelvin y: ppm *10^-6 / K
        piece_CTE.curvefit = @(x) (-7.57545e-18).*x.^4 + (4.34715e-14).*x.^3 + (-6.58705e-11).*x.^2 + (4.10566e-8).*x + (4.59464e-6);
        piece_CTE.extrapol = 16e-6;
        piece_CTE.extrastart = 1033.15;
        CTE = @(x) CTE_func(x,piece_CTE);

        % k curve x: Kelvin y: W/m K
     

        piece_k.ascending = @(T) 0.0143*T + 6.2118;
        piece_k.inflection = 1675.15;
        k_inflection = piece_k.ascending(piece_k.inflection);
        piece_k.decending = @(T) k_inflection - 0.0143.*(T-1675.15);
        k = @(T) k_func(T,piece_k);

        % Conductivity antiderivative
        piece_K_ad.ascending = @(T) (0.0143/2).*(T.^2) + 6.2118.*T;
        piece_K_ad.inflection = 1675.15;
        K_inflection = piece_K_ad.ascending(piece_K_ad.inflection);
        piece_K_ad.decending = @(T) K_inflection + k_inflection.*(T-1675.15) - (0.0143/2).*((T-1675.15).^2);
        K_ad = @(T) K_ad_func(T,piece_K_ad);

        % Youngs Modulus [Pa]
        piece_E.interp = @(T) interp1(E_data(:,1), E_data(:,2), T, 'linear');
        piece_E.lowertempbound = E_data(1,1);
        piece_E.uppertempbound = E_data(end,1);
        piece_E.lowerbound = E_data(1,2);
        piece_E.upperbound = E_data(end,2);
        E = @(x) E_func(x,piece_E);

        % Elongation [fractional] 
        piece_elong.interp = @(T) interp1(elong_data(:,1), elong_data(:,2), T, 'linear');
        piece_elong.lowertempbound = elong_data(1,1);
        piece_elong.uppertempbound = elong_data(end,1);
        piece_elong.lowerbound = elong_data(1,2);
        piece_elong.upperbound = elong_data(end,2);
        elong = @(x) elong_func(x,piece_elong);

        % Yield Strength [Mpa]
        piece_yts.interp = @(T) interp1(yts_data(:,1), yts_data(:,2), T, 'linear');
        piece_yts.lowertempbound = yts_data(1,1);
        piece_yts.uppertempbound = yts_data(end,1);
        piece_yts.lowerbound = yts_data(1,2);
        piece_yts.upperbound = yts_data(end,2);
        yts = @(x) yts_func(x,piece_yts);
    elseif material == "grcop42"
        properties = readmatrix(pwd + "/Material Data/GrCop42.xlsx");
        k = properties(13:end,1:2); % thermal conductivity [W/m-K] 
        E = [properties(1:13, 9) properties(1:13,10)];
        CTE = [properties(1:5,1) properties(1:5,3)]; % [ppm] 
        yts = properties(1:8,1:2);
        elong = [properties(1:8,1) properties(1:8,5)];
    elseif material == "inco-625"
        properties = readmatrix(pwd + "/Material Data/Inconel625.xlsx");
        k = properties(8:19,1:2); % thermal conductivity [W/m-K] 
        E = [properties(1:5, 1) properties(1:5,5)];
        CTE = [properties(22:end,1:2)]; % [ppm] 
        yts = properties(1:5,1:2);
        elong = [properties(1:5,1) properties(1:5,3)];
    elseif material == "c101"
        properties = readmatrix(pwd + "/Material Data/C101.xlsx");
        k = properties(13:21,1:2); % thermal conductivity [W/m-K] 
        E = [properties(1:5, 9:10)];
        CTE = [properties(1:7,1) properties(1:7,3)]; % [ppm] 
        yts = properties(1:7, 1:2);
        elong = [properties(1:7,1) properties(1:7,5)];
    elseif material == "alsi10mg"
        properties = readmatrix(pwd + "/Material Data/AlSi10Mg.xlsx");
        k = properties(13:15,1:2); % thermal conductivity [W/m-K] 
        E = [properties(1:6, 9:10)];
        CTE = [properties(1:4,1) properties(1:4,3)]; % [ppm] 
        yts = properties(1:4,1:2);
        elong = [properties(1:4,1) properties(1:4,5)];
    end

    % Contruct material obj
    mat.CTE = CTE;
    mat.k = k;
    mat.K_ad = K_ad;
    mat.E = E;
    mat.elong = elong;
    mat.yts = yts;

end

function k = k_func(x,p)
    k = nan(size(x));

    idx = x <= p.inflection;
    k(idx) = p.ascending(x(idx));

    idx = x > p.inflection;
    k(idx) = p.decending(x(idx));

end

function K_ad = K_ad_func(x,p)
    K_ad = nan(size(x));

    idx = x <= p.inflection;
    K_ad(idx) = p.ascending(x(idx));

    idx = x > p.inflection;
    K_ad(idx) = p.decending(x(idx));

end

function CTE = CTE_func(x,p)
    CTE = nan(size(x));

    idx = x <= p.extrastart;
    CTE(idx) = p.curvefit(x(idx));

    idx = x > p.extrastart;
    CTE(idx) = p.extrapol;

end

function E = E_func(x,p)
    E = nan(size(x));

    idx = x <= p.lowertempbound;
    E(idx) = p.lowerbound;

    idx = x > p.lowertempbound & x < p.uppertempbound;
    E(idx) = p.interp(x(idx));

    idx = x >= p.uppertempbound;
    E(idx) = p.upperbound;

end

function elong = elong_func(x,p)
    elong = nan(size(x));

    idx = x <= p.lowertempbound;
    elong(idx) = p.lowerbound;

    idx = x > p.lowertempbound & x < p.uppertempbound;
    elong(idx) = p.interp(x(idx));

    idx = x >= p.uppertempbound;
    elong(idx) = p.upperbound;

end

function yts = yts_func(x,p)
    yts = nan(size(x));

    idx = x <= p.lowertempbound;
    yts(idx) = p.lowerbound;

    idx = x > p.lowertempbound & x < p.uppertempbound;
    yts(idx) = p.interp(x(idx));

    idx = x >= p.uppertempbound;
    yts(idx) = p.upperbound;

end

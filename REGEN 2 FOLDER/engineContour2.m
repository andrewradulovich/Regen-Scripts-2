%% Engine Contour 2
% Only does conical nozzles for now

function geom = engineContour2(R_t, exp_ratio, con_ratio, conv_angle, conical_half_angle, L_star, R_fillets)
%{ 
Description: 

Inputs:
- 

Outputs: 
- 
%}
    %% Existing geometry calculations

    A_t = pi*R_t^2; % Throat Area
    A_c = A_t*con_ratio; % Chamber Area

    R_c = sqrt(A_c/pi); % Chamber radius
    R_e = sqrt(exp_ratio)*R_t; % Exit radius
 
    theta_i  = conical_half_angle*pi/180; % nozzle expansion (initial) angle [rad]
    conv_angle = conv_angle*pi/180; % convergence angle [rad]

    Rcf = R_fillets(1)*R_t;   % chamber fillet radius
    Rtf = R_fillets(2)*R_t;   % converging throat fillet
    Rdf = R_fillets(3)*R_t;   % diverging throat fillet

    %% Define exact throat geometry

    % Converging throat arc
    x_tf_start = -Rtf*sin(conv_angle); % Start of converging radius section
    x_tf_end   = 0; % End of converging radius section

    r_tf_center = R_t + Rtf; % Height of the center of radius

    r_tf = @(x) r_tf_center - sqrt(Rtf^2 - x.^2); % R points of nozzle contour

    dr_tf = @(x) x ./ sqrt(Rtf^2 - x.^2); % Derivative w/ respect to x

    d2r_tf = @(x) Rtf^2 ./ (Rtf^2 - x.^2).^(3/2); % Second Derivative w/ respect to x


    % Diverging throat arc
    x_df_start = 0; %Start of diverging radius section
    x_df_end   = Rdf*sin(theta_i); % end of diverging radius section

    r_df_center = R_t + Rdf; % Height of the center of radius

    r_df = @(x) r_df_center - sqrt(Rdf^2 - x.^2); % R points of nozzle contour

    dr_df = @(x) x ./ sqrt(Rdf^2 - x.^2); % Derivative w/ respect to x 

    d2r_df = @(x) Rdf^2 ./ (Rdf^2 - x.^2).^(3/2); % Second Derivative w/ respect to x


    %% Chamber converging fillet

    % First determine its x-center exactly as your old geometry does

    r_chamber_arc_end = R_c - Rcf + Rcf*cos(conv_angle);

    r_throat_arc_start = R_t + Rtf - Rtf*cos(conv_angle);

    x_cc = x_tf_start - (r_chamber_arc_end - r_throat_arc_start)/tan(conv_angle) + Rcf*cos(pi/2 + conv_angle);

    x_cf_start = x_cc;
    x_cf_end   = x_cc + Rcf*sin(conv_angle);

    r_cf_center = R_c - Rcf;

    r_cf = @(x) r_cf_center  + sqrt(Rcf^2 - (x-x_cc).^2);

    dr_cf = @(x) -(x-x_cc) ./ sqrt(Rcf^2 - (x-x_cc).^2);

    d2r_cf = @(x) -Rcf^2 ./ (Rcf^2 - (x-x_cc).^2).^(3/2);


    %% Linear converging section

    x_cl_start = x_cf_end;
    x_cl_end   = x_tf_start;

    r_cl_start = r_cf(x_cl_start);
    r_cl_end   = r_tf(x_cl_end);

    m_conv = (r_cl_end-r_cl_start)  / (x_cl_end-x_cl_start);

    r_cl = @(x) r_cl_start + m_conv*(x-x_cl_start);

    dr_cl = @(x) m_conv + zeros(size(x));

    d2r_cl = @(x) zeros(size(x));


    %% Chamber length

    vol_chamber = L_star*A_t;

    % For now this could be calculated by integral()

    vol_chamberRadius = pi * integral(@(x) r_cf(x).^2, x_cf_start, x_cf_end);
    vol_convLinear = pi * integral(@(x) r_cl(x).^2, x_cl_start, x_cl_end);
    vol_convThroatRadius = pi * integral(@(x) r_tf(x).^2, x_tf_start, x_tf_end);



    vol_converging = vol_chamberRadius + vol_convLinear + vol_convThroatRadius;

    vol_cylindrical = vol_chamber-vol_converging;
    if vol_cylindrical <= 0
        error('L* is too small for the specified converging geometry.');
    end
    L_c = vol_cylindrical/(pi*R_c^2);

    x_chamber_start = x_cf_start-L_c;
    x_chamber_end   = x_cf_start;

    r_chamber = @(x) R_c + zeros(size(x));
    dr_chamber = @(x) zeros(size(x));
    d2r_chamber = @(x) zeros(size(x));


    %% Conical nozzle

    % Preserve old geometry definition exactly
    L_nozzle = (R_e-R_t)/tan(pi/12);

    x_cone_start = x_df_end;
    x_cone_end   = x_df_end + L_nozzle;

    r_cone_start = r_df(x_cone_start);

    m_cone = (R_e-r_cone_start) ...
           / (x_cone_end-x_cone_start);

    r_cone = @(x) r_cone_start ...
        + m_cone*(x-x_cone_start);

    dr_cone = @(x) m_cone + zeros(size(x));

    d2r_cone = @(x) zeros(size(x));


    %% Package pieces

    pieces.x_chamber_start = x_chamber_start;
    pieces.x_cf_start      = x_cf_start;
    pieces.x_cf_end        = x_cf_end;
    pieces.x_tf_start      = x_tf_start;
    pieces.x_throat        = 0;
    pieces.x_df_end        = x_df_end;
    pieces.x_exit          = x_cone_end;

    pieces.r_chamber = r_chamber;
    pieces.r_cf      = r_cf;
    pieces.r_cl      = r_cl;
    pieces.r_tf      = r_tf;
    pieces.r_df      = r_df;
    pieces.r_cone    = r_cone;

    pieces.dr_chamber = dr_chamber;
    pieces.dr_cf      = dr_cf;
    pieces.dr_cl      = dr_cl;
    pieces.dr_tf      = dr_tf;
    pieces.dr_df      = dr_df;
    pieces.dr_cone    = dr_cone;

    pieces.d2r_chamber = d2r_chamber;
    pieces.d2r_cf = d2r_cf;
    pieces.d2r_cl = d2r_cl;
    pieces.d2r_tf = d2r_tf;
    pieces.d2r_df = d2r_df;
    pieces.d2r_cone = d2r_cone;

    pieces.k_chamber = 0;
    pieces.k_cf = 1/Rcf;
    pieces.k_cl = 0;
    pieces.k_tf = -1/Rtf;
    pieces.k_df = -1/Rdf;
    pieces.k_cone = 0;

    geom.r    = @(x) radiusAtX(x,pieces);
    geom.drdx = @(x) slopeAtX(x,pieces);
    geom.d2rdx = @(x) d2rAtX(x,pieces);
    geom.kappa = @(x) curvatureAtX(x,pieces);
    geom.dxds = @(x) ((1 + (geom.drdx(x)).^2).^(1/2)).^(-1);
    geom.aprime = @(x) -(geom.drdx(x).*geom.d2rdx(x)) ./ (1 + geom.drdx(x).^2).^(3/2); % d(dxds)/dx

    %% Useful bounds

    geom.bounds.chamber =  [x_chamber_start, x_cf_start];

    geom.bounds.converging =  [x_cf_start, 0];

    geom.bounds.diverging =  [0, x_cone_end];

    % I would also expose the individual geometric pieces
    geom.bounds.chamberFillet = [x_cf_start, x_cf_end];

    geom.bounds.convergingLine = [x_cf_end, x_tf_start];

    geom.bounds.convergingThroat = [x_tf_start, 0];

    geom.bounds.divergingThroat = [0, x_df_end];

    geom.bounds.nozzleLine = [x_df_end, x_cone_end];

    geom.R_t = R_t;
    geom.R_c = R_c;
    geom.R_e = R_e;

end

function r = radiusAtX(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x < p.x_cf_start;
    r(idx) = p.r_chamber(x(idx));

    idx = x >= p.x_cf_start & x < p.x_cf_end;
    r(idx) = p.r_cf(x(idx));

    idx = x >= p.x_cf_end & x < p.x_tf_start;
    r(idx) = p.r_cl(x(idx));

    idx = x >= p.x_tf_start & x <= p.x_throat;
    r(idx) = p.r_tf(x(idx));

    idx = x > p.x_throat & x <= p.x_df_end;
    r(idx) = p.r_df(x(idx));

    idx = x > p.x_df_end & x <= p.x_exit;
    r(idx) = p.r_cone(x(idx));

end

function drdx = slopeAtX(x,p)

    drdx = nan(size(x));

    idx = x >= p.x_chamber_start & x < p.x_cf_start;
    drdx(idx) = p.dr_chamber(x(idx));

    idx = x >= p.x_cf_start & x < p.x_cf_end;
    drdx(idx) = p.dr_cf(x(idx));

    idx = x >= p.x_cf_end & x < p.x_tf_start;
    drdx(idx) = p.dr_cl(x(idx));

    idx = x >= p.x_tf_start & x <= p.x_throat;
    drdx(idx) = p.dr_tf(x(idx));

    idx = x > p.x_throat & x <= p.x_df_end;
    drdx(idx) = p.dr_df(x(idx));

    idx = x > p.x_df_end & x <= p.x_exit;
    drdx(idx) = p.dr_cone(x(idx));

end

function drdx = d2rAtX(x,p)

    drdx = nan(size(x));

    idx = x >= p.x_chamber_start & x < p.x_cf_start;
    drdx(idx) = p.d2r_chamber(x(idx));

    idx = x >= p.x_cf_start & x < p.x_cf_end;
    drdx(idx) = p.d2r_cf(x(idx));

    idx = x >= p.x_cf_end & x < p.x_tf_start;
    drdx(idx) = p.d2r_cl(x(idx));

    idx = x >= p.x_tf_start & x <= p.x_throat;
    drdx(idx) = p.d2r_tf(x(idx));

    idx = x > p.x_throat & x <= p.x_df_end;
    drdx(idx) = p.d2r_df(x(idx));

    idx = x > p.x_df_end & x <= p.x_exit;
    drdx(idx) = p.d2r_cone(x(idx));

end

function kappa = curvatureAtX(x,p)
    
    kappa = nan(size(x));

    idx = x >= p.x_chamber_start & x < p.x_cf_start;
    kappa(idx) = p.k_chamber;

    idx = x >= p.x_cf_start & x < p.x_cf_end;
    kappa(idx) = p.k_cf;

    idx = x >= p.x_cf_end & x < p.x_tf_start;
    kappa(idx) = p.k_cl;

    idx = x >= p.x_tf_start & x <= p.x_throat;
    kappa(idx) = p.k_tf;

    idx = x > p.x_throat & x <= p.x_df_end;
    kappa(idx) = p.k_df;

    idx = x > p.x_df_end & x <= p.x_exit;
    kappa(idx) = p.k_cone;

end
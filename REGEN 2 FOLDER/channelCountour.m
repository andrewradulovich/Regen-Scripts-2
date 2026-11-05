%% Channel Contour

function geom = channelCountour(geom_enginecontour, h_c, w_c, t_w, float_bool, h_c_float, w_c_float, t_w_float)
%{ 
Description: Creates channel contour for a given channel height and mid radius width.
% Channel geometry inputs are taken and then linearly interpolate and
returned as a piecewise function. This allows for future implementation of
more complicated channel geometry in the future throught this function

Inputs:
Channel Geometry, each parameter should be a vector of four locations
[chamber_start, chamber_end, converging start, throat, diverging end]. 
Additionally, one optional float point where the user can define the
location and channel geometry can be defined
- bounds: output of engineContour2 geom.bounds
- h_c: Height of channels [im] at four sections
- w_c: Width of channels at mid radius [in] at four sections
- t_w: thickness of thin wall [in] at four sections
- Float_bool: 0: no extra point, 1: extra point in converging section


Outputs: 
- tw(x) [in]
- AR(x) [in]
- w_c(x) [in]
- h_c(x) [in]
- 
%}
    bounds = geom_enginecontour.bounds;
    r_fun = geom_enginecontour.r;

    if float_bool 
        % locate float point
        r_chamber = geom_enginecontour.r(bounds.chamber(1));
        A_chamber = pi()*r_chamber^2; 
        r_float = sqrt((AR_float * A_chamber)/ pi());
        x_float = fzero(@(x) r_fun - r_float,[0,[bounds.converging(1),bounds.converging(2)]]);
        
        % Calculate Slope & Intercepts

        % wall thickness
        m_tw_chamber = (t_w(2) - t_w(1)) / (bounds.chamber(2) - bounds.chamber(1)); 
        b_tw_chamber = t_w(1) - m_tw_chamber*bounds.chamber(1);

        m_tw_converging_upstr = (t_w_float - t_w(2)) / (x_float - bounds.converging(1));
        b_tw_converging_upstr = t_w(2) - m_tw_converging_upstr*bounds.converging(1);

        m_tw_converging_downstr = (t_w(3) - t_w_float) / (bounds.converging(2) - x_float);
        b_tw_converging_downstr = t_w_float - m_tw_converging_downstr*x_float;

        m_tw_diverging = (t_w(4) - t_w(3)) / (bounds.diverging(2)- bounds.diverging(1));
        b_tw_diverging = t_w(3) - m_tw_diverging*bounds.diverging(1);

        % channel height
        m_hc_chamber = (h_c(2) - h_c(1)) / (bounds.chamber(2) - bounds.chamber(1)); 
        b_hc_chamber = h_c(1) - m_hc_chamber*bounds.chamber(1);

        m_hc_converging_upstr = (h_c_float - h_c(2)) / (x_float - bounds.converging(1));
        b_hc_converging_upstr = h_c(2) - m_hc_converging_upstr*bounds.converging(1);

        m_hc_converging_downstr = (h_c(3) - h_c_float) / (bounds.converging(2) - x_float);
        b_hc_converging_downstr = h_c_float - m_hc_converging_downstr*x_float;

        m_hc_diverging = (h_c(4) - h_c(3)) / (bounds.diverging(2)- bounds.diverging(1));
        b_hc_diverging = h_c(3) - m_hc_diverging*bounds.diverging(1);
        % channel width
        m_wc_chamber = (w_c(2) - w_c(1)) / (bounds.chamber(2) - bounds.chamber(1)); 
        b_wc_chamber = w_c(1) - m_wc_chamber*bounds.chamber(1);

        m_wc_converging_upstr = (w_c_float - w_c(2)) / (x_float - bounds.converging(1));
        b_wc_converging_upstr = w_c(2) - m_wc_converging_upstr*bounds.converging(1);

        m_wc_converging_downstr = (w_c(3) - w_c_float) / (bounds.converging(2) - x_float);
        b_wc_converging_downstr = w_c_float - m_wc_converging_downstr*x_float;

        m_wc_diverging = (w_c(4) - w_c(3)) / (bounds.diverging(2)- bounds.diverging(1));
        b_wc_diverging = w_c(3) - m_wc_diverging*bounds.diverging(1);

        % Create Piecewise function
        tw_chamber = @(x) m_tw_chamber*x + b_tw_chamber;
        tw_converging_upstr = @(x) m_tw_converging_upstr*x + b_tw_converging_upstr;
        tw_converging_downstr = @(x) m_tw_converging_downstr*x + b_tw_converging_downstr;
        tw_diverging = @(x) m_tw_diverging*x + b_tw_diverging;

        dtwdx_chamber = m_tw_chamber;
        dtwdx_converging_upstr = m_tw_converging_upstr;
        dtwdx_converging_downstr = m_tw_converging_downstr;
        dtwdx_diverging = m_tw_diverging;

        hc_chamber = @(x) m_hc_chamber*x + b_hc_chamber;
        hc_converging_upstr = @(x) m_hc_converging_upstr*x + b_hc_converging_upstr;
        hc_converging_downstr = @(x) m_hc_converging_downstr*x + b_hc_converging_downstr;
        hc_diverging = @(x) m_hc_diverging*x + b_hc_diverging;

        dhcdx_chamber = m_hc_chamber;
        dhcdx_converging_upstr = m_hc_converging_upstr;
        dhcdx_converging_downstr = m_hc_converging_downstr;
        dhcdx_diverging = m_hc_diverging;

        wc_chamber = @(x) m_wc_chamber*x + b_wc_chamber;
        wc_converging_upstr = @(x) m_wc_converging_upstr*x + b_wc_converging_upstr;
        wc_converging_downstr = @(x) m_wc_converging_downstr*x + b_wc_converging_downstr;
        wc_diverging = @(x) m_wc_diverging*x + b_wc_diverging;

        dwcdx_chamber = m_wc_chamber;
        dwcdx_converging_upstr = m_wc_converging_upstr;
        dwcdx_converging_downstr = m_wc_converging_downstr;
        dwcdx_diverging = m_wc_diverging;

    
        % Put together pieces
        pieces.x_chamber_start = bounds.chamber(1);
        pieces.x_chamber_end = bounds.chamber(2);
        pieces.x_converging_upstr_start = bounds.converging(1);
        pieces.x_converging_upstr_end = x_float;
        pieces.x_converging_downstr_start = x_float;
        pieces.x_converging_downstr_end = bounds.converging(2);
        pieces.x_diverging_start = bounds.diverging(1);
        pieces.x_diverging_end = bounds.diverging(2); 

        pieces.tw_chamber = tw_chamber;
        pieces.tw_converging_upstr = tw_converging_upstr;
        pieces.tw_converging_downstr = tw_converging_downstr;
        pieces.tw_diverging = tw_diverging;

        pieces.dtwdx_chamber = dtwdx_chamber;
        pieces.dtwdx_converging_upstr = dtwdx_converging_upstr;
        pieces.dtwdx_converging_downstr = dtwdx_converging_downstr;
        pieces.dtwdx_diverging = dtwdx_diverging;

        pieces.hc_chamber = hc_chamber;
        pieces.hc_converging_upstr = hc_converging_upstr;
        pieces.hc_converging_downstr = hc_converging_downstr;
        pieces.hc_diverging = hc_diverging;

        pieces.dhcdx_chamber = dhcdx_chamber;
        pieces.dhcdx_converging_upstr = dhcdx_converging_upstr;
        pieces.dhcdx_converging_downstr = dhcdx_converging_downstr;
        pieces.dhcdx_diverging = dhcdx_diverging;

        pieces.wc_chamber = wc_chamber;
        pieces.wc_converging_upstr = wc_converging_upstr;
        pieces.wc_converging_downstr = wc_converging_downstr;
        pieces.wc_diverging = wc_diverging;

        pieces.dwcdx_chamber = dwcdx_chamber;
        pieces.dwcdx_converging_upstr = dwcdx_converging_upstr;
        pieces.dwcdx_converging_downstr = dwcdx_converging_downstr;
        pieces.dwcdx_diverging = dwcdx_diverging;
        
        geom.tw = @(x) wallthicknessAtX1(x,pieces);
        geom.dtwdx = @(x) ddxwallthicknessAtX1(x,pieces);
        geom.hc = @(x) heightAtX1(x,pieces);
        geom.dhxdx = @(x) ddxheightAtX1(x,pieces);
        geom.wc = @(x) widthAtX1(x,pieces);
        geom.wc = @(x) ddxwidthAtX1(x,pieces);

    else 
        % Calculate Slope & Intercepts

        % wall thickness
        m_tw_chamber = (t_w(2) - t_w(1)) / (bounds.chamber(2) - bounds.chamber(1)); 
        b_tw_chamber = t_w(1) - m_tw_chamber*bounds.chamber(1);

        m_tw_converging = (t_w(3) - t_w(2)) / (bounds.converging(2) - bounds.converging(1));
        b_tw_converging = t_w(2) - m_tw_converging*bounds.converging(1);
 
        m_tw_diverging = (t_w(4) - t_w(3)) / (bounds.diverging(2)- bounds.diverging(1));
        b_tw_diverging = t_w(3) - m_tw_diverging*bounds.diverging(1);

        % channel height
        m_hc_chamber = (h_c(2) - h_c(1)) / (bounds.chamber(2) - bounds.chamber(1)); 
        b_hc_chamber = h_c(1) - m_hc_chamber*bounds.chamber(1);

        m_hc_converging = (h_c(3) - h_c(2)) / (bounds.converging(2) - bounds.converging(1));
        b_hc_converging = h_c(2) - m_hc_converging*bounds.converging(1);

        m_hc_diverging = (h_c(4) - h_c(3)) / (bounds.diverging(2)- bounds.diverging(1));
        b_hc_diverging = h_c(3) - m_hc_diverging*bounds.diverging(1);
        % channel width
        m_wc_chamber = (w_c(2) - w_c(1)) / (bounds.chamber(2) - bounds.chamber(1)); 
        b_wc_chamber = w_c(1) - m_wc_chamber*bounds.chamber(1);

        m_wc_converging = (w_c(3) - w_c(2)) / (bounds.converging(2) - bounds.converging(1));
        b_wc_converging = w_c(2) - m_wc_converging*bounds.converging(1);

        m_wc_diverging = (w_c(4) - w_c(3)) / (bounds.diverging(2)- bounds.diverging(1));
        b_wc_diverging = w_c(3) - m_wc_diverging*bounds.diverging(1);

        % Create Piecewise function
        tw_chamber = @(x) m_tw_chamber*x + b_tw_chamber;
        tw_converging = @(x) m_tw_converging*x + b_tw_converging;
        tw_diverging = @(x) m_tw_diverging*x + b_tw_diverging;

        dtwdx_chamber = m_tw_chamber;
        dtwdx_converging = m_tw_converging;
        dtwdx_diverging = m_tw_diverging; 

        hc_chamber = @(x) m_hc_chamber*x + b_hc_chamber;
        hc_converging = @(x) m_hc_converging*x + b_hc_converging;
        hc_diverging = @(x) m_hc_diverging*x + b_hc_diverging;

        dhcdx_chamber = m_hc_chamber;
        dhcdx_converging = m_hc_converging;
        dhcdx_diverging = m_hc_diverging;

        wc_chamber = @(x) m_wc_chamber*x + b_wc_chamber;
        wc_converging = @(x) m_wc_converging*x + b_wc_converging;
        wc_diverging = @(x) m_wc_diverging*x + b_wc_diverging;

        dwcdx_chamber = m_wc_chamber;
        dwcdx_converging = m_wc_converging;
        dwcdx_diverging = m_wc_diverging;
    
        % Put together pieces
        pieces.x_chamber_start = bounds.chamber(1);
        pieces.x_chamber_end = bounds.chamber(2);
        pieces.x_converging_start = bounds.converging(1);
        pieces.x_converging_end = bounds.converging(2);
        pieces.x_diverging_start = bounds.diverging(1);
        pieces.x_diverging_end = bounds.diverging(2); 

        pieces.tw_chamber = tw_chamber;
        pieces.tw_converging = tw_converging;
        pieces.tw_diverging = tw_diverging;

        pieces.dtwdx_chamber = dtwdx_chamber;
        pieces.dtwdx_converging = dtwdx_converging;
        pieces.dtwdx_diverging = dtwdx_diverging;

        pieces.hc_chamber = hc_chamber;
        pieces.hc_converging = hc_converging;
        pieces.hc_diverging = hc_diverging;

        pieces.dhcdx_chamber = dhcdx_chamber;
        pieces.dhcdx_converging = dhcdx_converging;
        pieces.dhcdx_diverging = dhcdx_diverging;

        pieces.wc_chamber = wc_chamber;
        pieces.wc_converging = wc_converging;
        pieces.wc_diverging = wc_diverging;
        
        pieces.dwcdx_chamber = dwcdx_chamber;
        pieces.dwcdx_converging = dwcdx_converging;
        pieces.dwcdx_diverging = dwcdx_diverging;

        geom.tw = @(x) wallthicknessAtX2(x,pieces);
        geom.dtwdx = @(x) ddxwallthicknessAtX2(x,pieces);
        geom.hc = @(x) heightAtX2(x,pieces);
        geom.dhcdx = @(x) ddxheightAtX2(x,pieces);
        geom.wc = @(x) widthAtX2(x,pieces);
        geom.dwcdx = @(x) ddxwidthAtX2(x,pieces);
    end
    
end
    
    
    
    
function r = wallthicknessAtX1(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_converging_upstr_start;
    r(idx) = p.tw_chamber(x(idx));

    idx = x > p.x_converging_upstr_start & x < p.x_converging_upstr_end;
    r(idx) = p.tw_converging_upstr(x(idx));

    idx = x >= p.x_converging_downstr_start & x < p.x_converging_downstr_end;
    r(idx) = p.tw_converging_downstr(x(idx));

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.tw_diverging(x(idx));
end

function r = ddxwallthicknessAtX1(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_converging_upstr_start;
    r(idx) = p.dtwdx_chamber;

    idx = x > p.x_converging_upstr_start & x < p.x_converging_upstr_end;
    r(idx) = p.dtwdx_converging_upstr;

    idx = x >= p.x_converging_downstr_start & x < p.x_converging_downstr_end;
    r(idx) = p.dtwdx_converging_downstr;

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.dtwdx_diverging;
end

function r = heightAtX1(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.hc_chamber(x(idx));

    idx = x > p.x_converging_upstr_start & x < p.x_converging_upstr_end;
    r(idx) = p.hc_converging_upstr(x(idx));

    idx = x >= p.x_converging_downstr_start & x < p.x_converging_downstr_end;
    r(idx) = p.hc_converging_downstr(x(idx));

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.hc_diverging(x(idx));
end

function r = ddxheightAtX1(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.dhcdx_chamber;

    idx = x > p.x_converging_upstr_start & x < p.x_converging_upstr_end;
    r(idx) = p.dhcdx_converging_upstr;

    idx = x >= p.x_converging_downstr_start & x < p.x_converging_downstr_end;
    r(idx) = p.dhcdx_converging_downstr;

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.dhcdx_diverging;
end

function r = widthAtX1(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.wc_chamber(x(idx));

    idx = x > p.x_converging_upstr_start & x < p.x_converging_upstr_end;
    r(idx) = p.wc_converging_upstr(x(idx));

    idx = x >= p.x_converging_downstr_start & x < p.x_converging_downstr_end;
    r(idx) = p.wc_converging_downstr(x(idx));

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.wc_diverging(x(idx));
end

function r = ddxwidthAtX1(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.dwcdx_chamber;

    idx = x > p.x_converging_upstr_start & x < p.x_converging_upstr_end;
    r(idx) = p.dwcdx_converging_upstr;

    idx = x >= p.x_converging_downstr_start & x < p.x_converging_downstr_end;
    r(idx) = p.dwcdx_converging_downstr;

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.dwcdx_diverging;
end

function r = wallthicknessAtX2(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.tw_chamber(x(idx));

    idx = x > p.x_converging_start & x < p.x_converging_end;
    r(idx) = p.tw_converging(x(idx));

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.tw_diverging(x(idx));
end

function r = ddxwallthicknessAtX2(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.dtwdx_chamber;

    idx = x > p.x_converging_start & x < p.x_converging_end;
    r(idx) = p.dtwdx_converging;

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.dtwdx_diverging;
end

function r = heightAtX2(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.hc_chamber(x(idx));

    idx = x > p.x_converging_start & x < p.x_converging_end;
    r(idx) = p.hc_converging(x(idx));

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.hc_diverging(x(idx));
end

function r = ddxheightAtX2(x,p)
    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.dhcdx_chamber;

    idx = x > p.x_converging_start & x < p.x_converging_end;
    r(idx) = p.dhcdx_converging;

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.dhcdx_diverging;
end

function r = widthAtX2(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.wc_chamber(x(idx));

    idx = x > p.x_converging_start & x < p.x_converging_end;
    r(idx) = p.wc_converging(x(idx));

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.wc_diverging(x(idx));
end

function r = ddxwidthAtX2(x,p)

    r = nan(size(x));

    idx = x >= p.x_chamber_start & x <= p.x_chamber_end;
    r(idx) = p.dwcdx_chamber;

    idx = x > p.x_converging_start & x < p.x_converging_end;
    r(idx) = p.dwcdx_converging;

    idx = x >= p.x_diverging_start & x <= p.x_diverging_end;
    r(idx) = p.dwcdx_diverging;
end
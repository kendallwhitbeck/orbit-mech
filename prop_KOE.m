% Function that takes 1x6 vector of initial orbital elements
% (SMA,Ecc,Inc,RAAN,ArgPeri,TruAnom) & converts it to Mx7 matrix where each
% row corresponds to the next time value given in 't_vec'.
function [KOE_mat] = prop_KOE(sma, ecc, inc, raan, arg, t_p, t_vec, mu)

for ii = 1:length(t_vec)
    
    % Calculates Mean Anomaly for elliptical orbit
    Me = sqrt(mu/sma^3)*(t_vec(ii)-t_p);
    
    % Ensure that Mean Anomaly remains under 2*pi radians
    while Me > 2*pi
        Me = Me - 2*pi;
    end
    
    % Call function converting mean anomaly to true anomaly
    tru = Me2tru(ecc,Me);
    
    % Assigning KOEs to KOE_mat
    KOE_mat(ii,:) = [sma ecc inc raan arg tru Me];
    
end
end
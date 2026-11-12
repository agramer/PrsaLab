function d = concent_param(alpha, degFlag)

if degFlag
    alpha = circ_ang2rad(alpha);
end

alpha = alpha(:);   % ensure column vectors

%--- mean direction ---
mu = angle(mean(exp(1i * alpha)));

%--- absolute angular deviation from the group mean ---
d = abs(angle(exp(1i * (alpha - mu))));   % in [0, pi]

d=mean(d);


function [time] = getTimeScales(s, u)
% From: v = dz/dt --> dt = dz/v --> t = int(dt) = int (dz/v)
    time = trapz(s,1./u);
end
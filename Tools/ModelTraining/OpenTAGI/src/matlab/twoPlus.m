function [m, S] = twoPlus(m, S, deltaM, deltaS)
max_delta_m = sqrt(S)/20;
max_delta_S = max_delta_m/20;
m = m + sign(deltaM).*min(max_delta_m,1*abs(deltaM));
S = S + sign(deltaS).*min(max_delta_S,abs(deltaS));
% m = m + deltaM;
% S = S + deltaS;
end
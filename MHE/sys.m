%% Function: sys (To calculate one step ahed prediction)
function X_sys = sys(Xc)
global ts A B

X_sys =  A*Xc(:,1);
end

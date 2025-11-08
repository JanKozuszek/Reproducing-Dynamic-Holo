
# ==============================================
#  Auto-generated Julia code from Mathematica
#  Do not edit manually
# ==============================================


# ----------------------------------------------
#  Generated functions
# ----------------------------------------------

function DS0f(t)

res = E^((H*t)/2)*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om)); 

return res;
end

function DS1f(t)

res = (E^((H*t)/2)*H*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(1 + tanh(Om*(t - tstar))))/2; 

return res;
end

function DS2f(t)

res = (E^((H*t)/2)*H*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(1 + tanh(Om*(t - tstar)))*(H + 2*Om + (H - 2*Om)*tanh(Om*(t - tstar))))/4; 

return res;
end

function DS3f(t)

res = (E^((H*t)/2)*H*sech(Om*(t - tstar))^2*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(2*(3*H - 2*Om)*Om + (H^2 + 4*Om^2)*cosh(2*Om*(t - tstar)) + (H^2 - 4*Om^2)*sinh(2*Om*(t - tstar)))*(1 + tanh(Om*(t - tstar))))/8; 

return res;
end

function DS4f(t)

res = (E^((H*t)/2)*H*sech(Om*(t - tstar))^2*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(1 + tanh(Om*(t - tstar)))*(-H^3 + 12*H^2*Om + 12*H*Om^2 - 32*Om^3 + 2*(H^3 + 8*Om^3)*cosh(2*Om*(t - tstar)) + (H^3 - 8*Om^3)*sech(Om*(t - tstar))*sinh(3*Om*(t - tstar)) + 12*H^2*Om*tanh(Om*(t - tstar)) - 44*H*Om^2*tanh(Om*(t - tstar)) + 40*Om^3*tanh(Om*(t - tstar))))/16; 

return res;
end

function DS5f(t)

res = (E^((H*t)/2)*H*sech(Om*(t - tstar))^4*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(100*H^2*Om^2 - 240*H*Om^3 + 176*Om^4 + 4*Om*(5*H^3 - 10*H^2*Om + 40*H*Om^2 - 48*Om^3)*cosh(2*Om*(t - tstar)) + (H^4 + 16*Om^4)*cosh(4*Om*(t - tstar)) + 20*H^3*Om*sinh(2*Om*(t - tstar)) - 40*H^2*Om^2*sinh(2*Om*(t - tstar)) - 80*H*Om^3*sinh(2*Om*(t - tstar)) + 160*Om^4*sinh(2*Om*(t - tstar)) + H^4*sinh(4*Om*(t - tstar)) - 16*Om^4*sinh(4*Om*(t - tstar)))*(1 + tanh(Om*(t - tstar))))/32; 

return res;
end

function DS6f(t)

res = (E^((H*t)/2)*H*sech(Om*(t - tstar))^6*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(cosh(Om*(t - tstar)) + sinh(Om*(t - tstar)))*(20*Om^2*(13*H^3 - 12*H^2*Om - 44*H*Om^2 + 64*Om^3)*cosh(Om*(t - tstar)) + 10*Om*(3*H^4 - 8*H^3*Om + 12*H^2*Om^2 + 40*H*Om^3 - 80*Om^4)*cosh(3*Om*(t - tstar)) + H^5*cosh(5*Om*(t - tstar)) + 32*Om^5*cosh(5*Om*(t - tstar)) + 260*H^3*Om^2*sinh(Om*(t - tstar)) - 1680*H^2*Om^3*sinh(Om*(t - tstar)) + 3792*H*Om^4*sinh(Om*(t - tstar)) - 2944*Om^5*sinh(Om*(t - tstar)) + 30*H^4*Om*sinh(3*Om*(t - tstar)) - 80*H^3*Om^2*sinh(3*Om*(t - tstar)) + 120*H^2*Om^3*sinh(3*Om*(t - tstar)) - 592*H*Om^4*sinh(3*Om*(t - tstar)) + 864*Om^5*sinh(3*Om*(t - tstar)) + H^5*sinh(5*Om*(t - tstar)) - 32*Om^5*sinh(5*Om*(t - tstar))))/64; 

return res;
end


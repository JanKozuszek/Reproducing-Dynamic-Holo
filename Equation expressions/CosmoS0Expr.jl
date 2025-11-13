
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

function DS1f(params)
    t = params[1];

res = (E^((H*t)/2)*H*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(1 + tanh(Om*(t - tstar))))/2; 

return res;
end

function DS2f(params, DS1, p3, p4)
t = params[1];

res = (E^((H*t)/2)*H*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(1 + tanh(Om*(t - tstar)))*(H + 2*Om + (H - 2*Om)*tanh(Om*(t - tstar))))/4; 

return res;
end

function DS3f(params, DS1, p3, p4)
    t = params[1];

res = (E^((H*t)/2)*H*sech(Om*(t - tstar))^2*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(2*(3*H - 2*Om)*Om + (H^2 + 4*Om^2)*cosh(2*Om*(t - tstar)) + (H^2 - 4*Om^2)*sinh(2*Om*(t - tstar)))*(1 + tanh(Om*(t - tstar))))/8; 

return res;
end

function DS4f(params, DS1, p3, p4)
    t = params[1];

res = (E^((H*t)/2)*H*sech(Om*(t - tstar))^2*(cosh(Om*(t - tstar))*sech(Om*tstar))^(H/(2*Om))*(1 + tanh(Om*(t - tstar)))*(-H^3 + 12*H^2*Om + 12*H*Om^2 - 32*Om^3 + 2*(H^3 + 8*Om^3)*cosh(2*Om*(t - tstar)) + (H^3 - 8*Om^3)*sech(Om*(t - tstar))*sinh(3*Om*(t - tstar)) + 12*H^2*Om*tanh(Om*(t - tstar)) - 44*H*Om^2*tanh(Om*(t - tstar)) + 40*Om^3*tanh(Om*(t - tstar))))/16; 

return res;
end

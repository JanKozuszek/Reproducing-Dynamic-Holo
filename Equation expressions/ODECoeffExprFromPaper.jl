
# ==============================================
#  Auto-generated Julia code from Mathematica
#  Do not edit manually
# ==============================================


# ----------------------------------------------
#  Generated functions
# ----------------------------------------------

function SCoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 12;
 else 
    co = (648 + (9*(6*DS0^3*Phi*z^3 + 2*DS0^3*PhiZ*z^4 + z*(-2*DS1^3*(1 - 4*LN)*z^3 + DS0^3*(2 - 4*X*z) + DS0*DS1*z^2*(DS2*(1 - 4*LN)*z + DS1*(1 - 3*LN - 3*X*(z - 4*LN*z))) + DS0^2*z^2*(DS3*(1 - 4*LN)*z + DS2*(1 - 3*LN - 3*X*(z - 4*LN*z)))))^2)/DS0^6)/54;
end


    return co
end

function SCoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 8*z;
end


    return co
end

function SCoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = z^2;
end


    return co
end

function SSrc(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = (-4*DS0^2 - 9*DS1^2 - 9*DS0*DS2 + 72*DS0^2*Phi - 48*DS0*DS1*X)/(18*DS0);
 else 
    co = ((2*z^4*(1 - 2*X*z + 3*Phi*z^2 + PhiZ*z^3 - LN*(3*(DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^2 + (2*(-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^3)/(3*DS0^3)) + ((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3))/z)^2*(DS1 + DS0*X + DS0/z - (DS0*z)/3 + (-1/9*DS1 + (DS0*X)/3)*z^2 - LN*((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5)))/3 + z^4*(2*(-1/9*DS1 + (DS0*X)/3) + (2*DS0)/z^3 - LN*(6*(-1/6*DS1^2/DS0 - DS2/6)*z + (2*(5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^2)/(5*DS0^2) - (9*(DS1^2 + DS0*DS2)^2*z^3)/(10*DS0^3) + 20*(((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^3) + (2*(3*(-1/6*DS1^2/DS0 - DS2/6)*z^2 + (2*(5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^3)/(15*DS0^2) - ((DS1^2 + DS0*DS2)^2*z^4)/(10*DS0^3) + 5*(((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^4))/z - ((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5)/z^2) + 2*z^3*(-1/3*DS0 - DS0/z^2 + 2*(-1/9*DS1 + (DS0*X)/3)*z - LN*(3*(-1/6*DS1^2/DS0 - DS2/6)*z^2 + (2*(5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^3)/(15*DS0^2) - ((DS1^2 + DS0*DS2)^2*z^4)/(10*DS0^3) + 5*(((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^4) + ((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5)/z))/z^5;
end


    return co
end

function SdotCoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 1;
 else 
    co = -8*S*z^4 + ((23*DS1^2 + 28*DS0*DS2)*(1 - 6*LN)*z^6)/(135*DS0) + (2*(DS1^2 + DS0*DS2)*(1 - 6*LN)*p2*z^6)/(5*DS0) - 2*z*(DS1 + DS0*X + SZ*z^4) + (z^2*(DS1^4*(199 - 1554*LN + 1080*LN^2)*z^4 - 600*DS0^4*(-2 + 3*X*z) + 3*DS0*DS1^2*z^3*(3*DS2*(-53 + 238*LN + 240*LN^2)*z + 100*DS1*(-1 + 5*LN + 4*X*(z - 6*LN*z))) + DS0^2*z^2*(2*DS2^2*(247 - 1662*LN + 540*LN^2)*z^2 + 60*DS1^2*(5 - 20*LN + 24*(1 - 6*LN)*X^2*z^2 - 15*X*(z - 5*LN*z)) - 5*DS1*z*(-47*DS3*(1 - 6*LN)*z + 84*DS2*(-1 + 5*LN + 4*X*(z - 6*LN*z)))) + 5*DS0^3*z*(120*DS1 + z*(12*DS2*(5 - 20*LN + 24*(1 - 6*LN)*X^2*z^2 - 15*X*(z - 5*LN*z)) + z*(25*DS4*(1 - 6*LN)*z - 48*DS3*(-1 + 5*LN + 4*X*(z - 6*LN*z)))))))/(900*DS0^3);
end


    return co
end

function SdotCoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = (z*(-5400*DS0^3*S*z^4 - 3*DS1^4*(199 - 180*LN)*LN*z^6 - 1800*DS0^4*(3 + 3*X*z + z^2*(-1 + X*z)) - 9*DS0*DS1^2*LN*z^5*(3*DS2*(-53 - 40*LN)*z + 100*DS1*(-1 + 4*X*z)) - DS0^2*LN*z^4*(460*DS1^2*z^2 + 1080*DS1^2*p2*z^2 + 3*(2*DS2^2*(247 - 90*LN)*z^2 + 60*DS1^2*(5 - 15*X*z + 24*X^2*z^2) - 5*DS1*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z)))) + 5*DS0^3*(120*DS1*z*(-9 + z^2) - LN*z^4*(112*DS2*z^2 + 216*DS2*p2*z^2 + 3*(12*DS2*(5 - 15*X*z + 24*X^2*z^2) + z*(25*DS4*z - 48*DS3*(-1 + 4*X*z)))))))/(5400*DS0^3);
end


    return co
end

function SdotCoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 0;
end


    return co
end

function SdotSrc(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = (20*DS0^2 - 72*a4*DS0^2 - 27*DS1^2 + 61*DS0*DS2 - 72*DS0^2*p2 - 32*DS0*DS1*X + 72*DS0^2*X^2)/(144*DS0);
 else 
    co = (-(z^2*((-2*DS1)/9 - DS0/z^3 - (DS1 + DS0*X)/z^2 - LN*(((36*DS1^2 - 12*DS0*DS2)*z)/(72*DS0) + ((-900*DS1^3 + 1620*DS0*DS1*DS2 - 360*DS0^2*DS3 - 2700*DS0*DS1^2*X + 900*DS0^2*DS2*X)*z^2)/(1800*DS0^2)) + (((36*DS1^2 - 12*DS0*DS2)*z^2)/(144*DS0) + ((-900*DS1^3 + 1620*DS0*DS1*DS2 - 360*DS0^2*DS3 - 2700*DS0*DS1^2*X + 900*DS0^2*DS2*X)*z^3)/(5400*DS0^2))/z)*(DS1 + DS0*X + DS0/z - (DS0*z)/3 + (-1/9*DS1 + (DS0*X)/3)*z^2 + S*z^3 - LN*((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5))) - 2*z^2*(-1/6*DS0 + DS1^2/(2*DS0) + DS1*X + (DS0*X^2)/2 + DS0/(2*z^2) + (DS1 + DS0*X)/z - (2*DS1*z)/9 - LN*(((36*DS1^2 - 12*DS0*DS2)*z^2)/(144*DS0) + ((-900*DS1^3 + 1620*DS0*DS1*DS2 - 360*DS0^2*DS3 - 2700*DS0*DS1^2*X + 900*DS0^2*DS2*X)*z^3)/(5400*DS0^2)))*(-1/3*DS0 - DS0/z^2 + 2*(-1/9*DS1 + (DS0*X)/3)*z + 3*S*z^2 + SZ*z^3 - LN*(3*(-1/6*DS1^2/DS0 - DS2/6)*z^2 + (2*(5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^3)/(15*DS0^2) - ((DS1^2 + DS0*DS2)^2*z^4)/(10*DS0^3) + 5*(((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^4) + ((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5)/z) + (2*(DS1 + DS0*X + DS0/z - (DS0*z)/3 + (-1/9*DS1 + (DS0*X)/3)*z^2 + S*z^3 - LN*((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5))^2*((-z + X*z^2 - Phi*z^3 + LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)) + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^3/4)^2/2 - (4*(-3/2 - (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^2/2 + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^4/16)^2)/3))/3)/z^2;
end


    return co
end

function PhidotCoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = -1/3*DS0;
 else 
    co = (-70200*DS0^3*S*z^4 + 3*DS1^4*(597 - 4861*LN + 3420*LN^2)*z^6 - 1800*DS0^4*(3 + 12*X*z + z^2*(-7 + 10*X*z)) + 9*DS0*DS1^2*z^5*(3*DS2*(-159 + 767*LN + 760*LN^2)*z + 100*DS1*(-3 + 16*LN + 4*(3 - 19*LN)*X*z)) + DS0^2*z^4*(460*DS1^2*(3 - 19*LN)*z^2 + 1080*DS1^2*(3 - 19*LN)*p2*z^2 + 3*(2*DS2^2*(741 - 5233*LN + 1710*LN^2)*z^2 + 60*DS1^2*(15 - 65*LN - 15*(3 - 16*LN)*X*z + 24*(3 - 19*LN)*X^2*z^2) + 5*DS1*z*(47*DS3*(3 - 19*LN)*z - 84*DS2*(-3 + 16*LN + 4*(3 - 19*LN)*X*z)))) + 5*DS0^3*(240*DS1*z*(-18 + 5*z^2) + z^4*(-3240*SZ*z + 112*DS2*(3 - 19*LN)*z^2 + 216*DS2*(3 - 19*LN)*p2*z^2 + 3*(12*DS2*(15 - 65*LN - 15*(3 - 16*LN)*X*z + 24*(3 - 19*LN)*X^2*z^2) + z*(25*DS4*(3 - 19*LN)*z - 48*DS3*(-3 + 16*LN + 4*(3 - 19*LN)*X*z))))))/(16200*DS0^3);
end


    return co
end

function PhidotCoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = (z*(-5400*DS0^3*S*z^4 - 3*DS1^4*(199 - 180*LN)*LN*z^6 - 1800*DS0^4*(3 + 3*X*z + z^2*(-1 + X*z)) - 9*DS0*DS1^2*LN*z^5*(3*DS2*(-53 - 40*LN)*z + 100*DS1*(-1 + 4*X*z)) - DS0^2*LN*z^4*(460*DS1^2*z^2 + 1080*DS1^2*p2*z^2 + 3*(2*DS2^2*(247 - 90*LN)*z^2 + 60*DS1^2*(5 - 15*X*z + 24*X^2*z^2) - 5*DS1*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z)))) + 5*DS0^3*(120*DS1*z*(-9 + z^2) - LN*z^4*(112*DS2*z^2 + 216*DS2*p2*z^2 + 3*(12*DS2*(5 - 15*X*z + 24*X^2*z^2) + z*(25*DS4*z - 48*DS3*(-1 + 4*X*z)))))))/(8100*DS0^3);
end


    return co
end

function PhidotCoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 0;
end


    return co
end

function PhidotSrc(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = DS0/9 - DS1^2/(4*DS0) + DS2/4 - (DS0*Phi)/2 + (DS0*X^2)/2;
 else 
    co = (-(z^2*(-1/6*DS0 + DS1^2/(2*DS0) + DS1*X + (DS0*X^2)/2 + DS0/(2*z^2) + (DS1 + DS0*X)/z - (2*DS1*z)/9 + Sdot*z^2 - LN*(((36*DS1^2 - 12*DS0*DS2)*z^2)/(144*DS0) + ((-900*DS1^3 + 1620*DS0*DS1*DS2 - 360*DS0^2*DS3 - 2700*DS0*DS1^2*X + 900*DS0^2*DS2*X)*z^3)/(5400*DS0^2)))*(1 - 2*X*z + 3*Phi*z^2 + PhiZ*z^3 - LN*(3*(DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^2 + (2*(-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^3)/(3*DS0^3)) + ((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3))/z)) - z^2*(-1/2 - LN*(((-9*DS1^2 - 9*DS0*DS2)*z^2)/(12*DS0^2) + ((12*DS1^3 - 6*DS0*DS1*DS2 - 6*DS0^2*DS3 + 18*DS0*DS1^2*X + 18*DS0^2*DS2*X)*z^3)/(12*DS0^3)))*(-1/3*DS0 - DS0/z^2 + 2*(-1/9*DS1 + (DS0*X)/3)*z + 3*S*z^2 + SZ*z^3 - LN*(3*(-1/6*DS1^2/DS0 - DS2/6)*z^2 + (2*(5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^3)/(15*DS0^2) - ((DS1^2 + DS0*DS2)^2*z^4)/(10*DS0^3) + 5*(((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^4) + ((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5)/z) - ((DS1 + DS0*X + DS0/z - (DS0*z)/3 + (-1/9*DS1 + (DS0*X)/3)*z^2 + S*z^3 - LN*((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5))*(2*z^2*(-(LN*(((-9*DS1^2 - 9*DS0*DS2)*z)/(6*DS0^2) + ((12*DS1^3 - 6*DS0*DS1*DS2 - 6*DS0^2*DS3 + 18*DS0*DS1^2*X + 18*DS0^2*DS2*X)*z^2)/(4*DS0^3))) + (((-9*DS1^2 - 9*DS0*DS2)*z^2)/(12*DS0^2) + ((12*DS1^3 - 6*DS0*DS1*DS2 - 6*DS0^2*DS3 + 18*DS0*DS1^2*X + 18*DS0^2*DS2*X)*z^3)/(12*DS0^3))/z) + (-1 + (3*(z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^2)/4)*(-z + X*z^2 - Phi*z^3 + LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)) + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^3/4) - (8*(-z + X*z^2 - Phi*z^3 + LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)) + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^3/4)*(-3/2 - (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^2/2 + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^4/16))/3))/3)/z^2;
end


    return co
end

function ACoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 6*DS0^2;
 else 
    co = (-5400*DS0^3*S*z^4 - 3*DS1^4*(199 - 180*LN)*LN*z^6 - 1800*DS0^4*(3 + 3*X*z + z^2*(-1 + X*z)) - 9*DS0*DS1^2*LN*z^5*(3*DS2*(-53 - 40*LN)*z + 100*DS1*(-1 + 4*X*z)) - DS0^2*LN*z^4*(460*DS1^2*z^2 + 1080*DS1^2*p2*z^2 + 3*(2*DS2^2*(247 - 90*LN)*z^2 + 60*DS1^2*(5 - 15*X*z + 24*X^2*z^2) - 5*DS1*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z)))) + 5*DS0^3*(120*DS1*z*(-9 + z^2) - LN*z^4*(112*DS2*z^2 + 216*DS2*p2*z^2 + 3*(12*DS2*(5 - 15*X*z + 24*X^2*z^2) + z*(25*DS4*z - 48*DS3*(-1 + 4*X*z))))))^2/(4860000*DS0^6);
end


    return co
end

function ACoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = (z*(-5400*DS0^3*S*z^4 - 3*DS1^4*(199 - 180*LN)*LN*z^6 - 1800*DS0^4*(3 + 3*X*z + z^2*(-1 + X*z)) - 9*DS0*DS1^2*LN*z^5*(3*DS2*(-53 - 40*LN)*z + 100*DS1*(-1 + 4*X*z)) - DS0^2*LN*z^4*(460*DS1^2*z^2 + 1080*DS1^2*p2*z^2 + 3*(2*DS2^2*(247 - 90*LN)*z^2 + 60*DS1^2*(5 - 15*X*z + 24*X^2*z^2) - 5*DS1*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z)))) + 5*DS0^3*(120*DS1*z*(-9 + z^2) - LN*z^4*(112*DS2*z^2 + 216*DS2*p2*z^2 + 3*(12*DS2*(5 - 15*X*z + 24*X^2*z^2) + z*(25*DS4*z - 48*DS3*(-1 + 4*X*z))))))^2)/(4860000*DS0^6);
end


    return co
end

function ACoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = (z^2*(-5400*DS0^3*S*z^4 - 3*DS1^4*(199 - 180*LN)*LN*z^6 - 1800*DS0^4*(3 + 3*X*z + z^2*(-1 + X*z)) - 9*DS0*DS1^2*LN*z^5*(3*DS2*(-53 - 40*LN)*z + 100*DS1*(-1 + 4*X*z)) - DS0^2*LN*z^4*(460*DS1^2*z^2 + 1080*DS1^2*p2*z^2 + 3*(2*DS2^2*(247 - 90*LN)*z^2 + 60*DS1^2*(5 - 15*X*z + 24*X^2*z^2) - 5*DS1*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z)))) + 5*DS0^3*(120*DS1*z*(-9 + z^2) - LN*z^4*(112*DS2*z^2 + 216*DS2*p2*z^2 + 3*(12*DS2*(5 - 15*X*z + 24*X^2*z^2) + z*(25*DS4*z - 48*DS3*(-1 + 4*X*z))))))^2)/(29160000*DS0^6);
end


    return co
end

function ASrc(v::PDEVars{T}) where {T}

    @unpack Phi, PhiZ, PhiZZ, S, SZ, SZZ, Sdot, SdotZ, SdotZZ, Phidot, PhidotZ, PhidotZZ, A, AZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = (-2*(4*DS0^2 + 9*DS1^2 + 15*DS0*DS2 - 45*DS0^2*Phi + 18*DS0^2*Phidot - 135*DS0*S + 54*DS0*Sdot + 18*DS0*DS1*X))/9;
 else 
    co = (12*z^2*(-1/6*DS0 + DS1^2/(2*DS0) + DS1*X + (DS0*X^2)/2 + DS0/(2*z^2) + (DS1 + DS0*X)/z - (2*DS1*z)/9 + Sdot*z^2 - LN*(((36*DS1^2 - 12*DS0*DS2)*z^2)/(144*DS0) + ((-900*DS1^3 + 1620*DS0*DS1*DS2 - 360*DS0^2*DS3 - 2700*DS0*DS1^2*X + 900*DS0^2*DS2*X)*z^3)/(5400*DS0^2)))*(-1/3*DS0 - DS0/z^2 + 2*(-1/9*DS1 + (DS0*X)/3)*z + 3*S*z^2 + SZ*z^3 - LN*(3*(-1/6*DS1^2/DS0 - DS2/6)*z^2 + (2*(5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^3)/(15*DS0^2) - ((DS1^2 + DS0*DS2)^2*z^4)/(10*DS0^3) + 5*(((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^4) + ((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5)/z) + (DS1 + DS0*X + DS0/z - (DS0*z)/3 + (-1/9*DS1 + (DS0*X)/3)*z^2 + S*z^3 - LN*((-1/6*DS1^2/DS0 - DS2/6)*z^3 + ((5*DS1^3 - 7*DS0*DS1*DS2 - 4*DS0^2*DS3 + 15*DS0*(DS1^2 + DS0*DS2)*X)*z^4)/(30*DS0^2) + (((DS1^2 + DS0*DS2)^2*LN)/(10*DS0^3) - (20*DS0^2*(23*DS1^2 + 28*DS0*DS2) + 1080*DS0^2*(DS1^2 + DS0*DS2)*p2 + 3*(199*DS1^4 - 477*DS0*DS1^2*DS2 + 494*DS0^2*DS2^2 + 235*DS0^2*DS1*DS3 + 125*DS0^3*DS4 - 240*DS0*(-5*DS1^3 + 7*DS0*DS1*DS2 + 4*DS0^2*DS3)*X + 1440*DS0^2*(DS1^2 + DS0*DS2)*X^2))/(5400*DS0^3))*z^5))^2*(z^4*(6/z^4 + (4*X)/z^3 - LN*((-4*DS2)/(3*DS0) + ((36*DS1*DS2 - 36*DS0*DS3 + 144*DS0*DS2*X)*z)/(18*DS0^2)) + (2*((-4*DS2*z)/(3*DS0) + ((36*DS1*DS2 - 36*DS0*DS3 + 144*DS0*DS2*X)*z^2)/(36*DS0^2)))/z - ((-2*DS2*z^2)/(3*DS0) + ((36*DS1*DS2 - 36*DS0*DS3 + 144*DS0*DS2*X)*z^3)/(108*DS0^2))/z^2) + 2*z^3*(-2/z^3 - (2*X)/z^2 - LN*((-4*DS2*z)/(3*DS0) + ((36*DS1*DS2 - 36*DS0*DS3 + 144*DS0*DS2*X)*z^2)/(36*DS0^2)) + ((-2*DS2*z^2)/(3*DS0) + ((36*DS1*DS2 - 36*DS0*DS3 + 144*DS0*DS2*X)*z^3)/(108*DS0^2))/z) - 4*z^2*(-1/2 + Phidot*z^2 - LN*(((-9*DS1^2 - 9*DS0*DS2)*z^2)/(12*DS0^2) + ((12*DS1^3 - 6*DS0*DS1*DS2 - 6*DS0^2*DS3 + 18*DS0*DS1^2*X + 18*DS0^2*DS2*X)*z^3)/(12*DS0^3)))*(1 - 2*X*z + 3*Phi*z^2 + PhiZ*z^3 - LN*(3*(DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^2 + (2*(-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^3)/(3*DS0^3)) + ((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3))/z) - (4*((-z + X*z^2 - Phi*z^3 + LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)) + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^3/4)^2/2 - (4*(-3/2 - (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^2/2 + (z - X*z^2 + Phi*z^3 - LN*((DS1^2/(2*DS0^2) + DS2/(2*DS0))*z^3 + ((-6*DS1^3 + 3*DS0*DS1*DS2 + 3*DS0^2*DS3 - 9*DS0*DS1^2*X - 9*DS0^2*DS2*X)*z^4)/(6*DS0^3)))^4/16)^2)/3))/3))/z^2;
end


    return co
end



# ==============================================
#  Auto-generated Julia code from Mathematica
#  Do not edit manually
# ==============================================


# ----------------------------------------------
#  Generated functions
# ----------------------------------------------

function SCoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 12;
 else 
    co = 12 + (2*DS1^3*(-1 + 4*LN)*M*z^4 + 2*DS0^3*z*(M - 2*M*X*z + 3*Phi*z^2 + PhiZ*z^3) + DS0*DS1*M*z^3*(DS2*(1 - 4*LN)*z + DS1*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)) + DS0^2*M*z^3*(DS3*(1 - 4*LN)*z + DS2*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)))^2/(6*DS0^6);
end


    return co
end

function SCoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 8*z;
end


    return co
end

function SCoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = z^2;
end


    return co
end

function SSrc(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = -1/18*(M*(9*DS1^2*M + 4*DS0^2*(M^3 - 18*p2) + DS0*(9*DS2*M + 48*DS1*M*X)))/DS0;
 else 
    co = ((2*DS1^3*(-1 + 4*LN)*M*z^3 + 2*DS0^3*(M - 2*M*X*z + 3*Phi*z^2 + PhiZ*z^3) + DS0*DS1*M*z^2*(DS2*(1 - 4*LN)*z + DS1*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)) + DS0^2*M*z^2*(DS3*(1 - 4*LN)*z + DS2*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)))^2*(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(-120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z))))) + 6*DS0^6*M*(-3*DS1^4*(2369 - 7950*LN + 2700*LN^2)*M*z^4 + 3600*DS0^4*M*(-1 + 3*X*z) + 9*DS0*DS1^2*M*z^3*(-3*DS2*(-543 + 1150*LN + 600*LN^2)*z + 100*DS1*(9 - 20*LN + 4*(-11 + 30*LN)*X*z)) + DS0^2*z^2*(-6*DS2^2*(2807 - 8400*LN + 1350*LN^2)*M*z^2 - 15*DS1*M*z*(47*DS3*(11 - 30*LN)*z + 84*DS2*(9 - 20*LN + 4*(-11 + 30*LN)*X*z)) + 20*DS1^2*(-135*(-9 + 20*LN)*M*X*z + 54*(-11 + 30*LN)*p2*z^2 + 216*(-11 + 30*LN)*M*X^2*z^2 + M*(-315 - 253*M^2*z^2 + 30*LN*(18 + 23*M^2*z^2)))) + 5*DS0^3*z*(-720*DS1*M + z*(3*M*z*(25*DS4*(-11 + 30*LN)*z - 48*DS3*(9 - 20*LN + 4*(-11 + 30*LN)*X*z)) + 4*DS2*(-135*(-9 + 20*LN)*M*X*z + 54*(-11 + 30*LN)*p2*z^2 + 216*(-11 + 30*LN)*M*X^2*z^2 + 60*LN*M*(9 + 14*M^2*z^2) - 7*M*(45 + 44*M^2*z^2))))))/(32400*DS0^9*z^2);
end


    return co
end

function SdotCoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 1;
 else 
    co = (2*(-3*DS1^4*(199 - 1374*LN + 540*LN^2)*M^2*z^5 + 9*DS0*DS1^2*M^2*z^4*(-3*DS2*(-53 + 278*LN + 120*LN^2)*z + 100*DS1*(1 - 5*LN + 4*(-1 + 6*LN)*X*z)) + 1800*DS0^4*(-2*M^2*z + 3*X*(1 + M^2*z^2)) + DS0^2*M*z^3*(-6*DS2^2*(247 - 1572*LN + 270*LN^2)*M*z^2 - 15*DS1*M*z*(47*DS3*(1 - 6*LN)*z + 84*DS2*(1 - 5*LN + 4*(-1 + 6*LN)*X*z)) + 20*DS1^2*(-135*(-1 + 5*LN)*M*X*z + 54*(-1 + 6*LN)*p2*z^2 + 216*(-1 + 6*LN)*M*X^2*z^2 + M*(-45 - 23*M^2*z^2 + 6*LN*(30 + 23*M^2*z^2)))) + 5*DS0^3*(4320*S*z^3 - 360*DS1*(-3 + M^2*z^2) + z^3*(1080*SZ*z + M*(3*M*z*(25*DS4*(-1 + 6*LN)*z - 48*DS3*(1 - 5*LN + 4*(-1 + 6*LN)*X*z)) + 4*DS2*(-135*(-1 + 5*LN)*M*X*z + 54*(-1 + 6*LN)*p2*z^2 + 216*(-1 + 6*LN)*M*X^2*z^2 + M*(-45 - 28*M^2*z^2 + 12*LN*(15 + 14*M^2*z^2))))))))/(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(1080*S*z^4 - 120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z)))));
end


    return co
end

function SdotCoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 1;
end


    return co
end

function SdotCoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 0;
end


    return co
end

function SdotSrc(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = -1/2*(a4*DS0) - (3*DS1^2*M^2)/(16*DS0) + (61*DS2*M^2)/144 + (5*DS0*M^4)/36 - (DS0*M*p2)/2 - (2*DS1*M^2*X)/9 + (DS0*M^2*X^2)/2;
 else 
    co = -1/8100*(45*DS0*(30*DS1^3*(1 - 3*LN)*M^2*z^5 + 180*DS0^3*(1 + X*z) - 9*DS0*DS1*M^2*z^4*(6*DS2*(1 - 3*LN)*z + 5*DS1*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)) + DS0^2*(20*DS1*z*(9 + 2*M^2*z^2) + 3*M^2*z^4*(4*DS3*(1 - 3*LN)*z + 5*DS2*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)))) + (90*DS0*(-30*DS1^3*LN*M^2*z^5 + 30*DS0^3*(-3 - 6*X*z + M^2*z^2 - 3*X^2*z^2) - 9*DS0*DS1*z^2*(-6*DS2*LN*M^2*z^3 + 5*DS1*(2 - LN*M^2*z^2 + 2*LN*M^2*X*z^3)) + DS0^2*(20*DS1*z*(-9 - 9*X*z + 2*M^2*z^2) + 3*LN*M^2*z^4*(-4*DS3*z + 5*DS2*(-1 + 2*X*z))))*(-3*DS1^4*(199 - 1175*LN + 450*LN^2)*M^2*z^6 + 1800*DS0^4*(-3 - M^2*z^2 + 2*M^2*X*z^3) + 9*DS0*DS1^2*M^2*z^5*(-3*DS2*(-53 + 225*LN + 100*LN^2)*z + 100*DS1*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + DS0^2*M*z^4*(-6*DS2^2*(247 - 1325*LN + 225*LN^2)*M*z^2 - 15*DS1*M*z*(47*DS3*(1 - 5*LN)*z + 84*DS2*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + 20*DS1^2*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*z^2 + 216*(-1 + 5*LN)*M*X^2*z^2 + M*(-45 - 23*M^2*z^2 + 5*LN*(27 + 23*M^2*z^2)))) + 5*DS0^3*z^3*(-240*DS1*M^2 + 3240*S*z + z*(1080*SZ*z + M*(3*M*z*(25*DS4*(-1 + 5*LN)*z - 48*DS3*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + 4*DS2*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*z^2 + 216*(-1 + 5*LN)*M*X^2*z^2 + M*(-45 - 28*M^2*z^2 + 5*LN*(27 + 28*M^2*z^2))))))))/(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(1080*S*z^4 - 120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z))))) + z^6*(3*DS1^4*(199 - 90*LN)*LN*M^2 - 9*DS0*DS1^2*LN*M^2*(3*DS2*(53 + 20*LN) + 100*DS1*(-4*X + z^(-1))) + DS0^2*LN*M*(6*DS2^2*(247 - 45*LN)*M + 15*DS1*M*(47*DS3 + 84*DS2*(-4*X + z^(-1))) + 20*DS1^2*(23*M^3 + 54*p2 + 216*M*X^2 + (45*M)/z^2 - (135*M*X)/z)) + 5*DS0^3*(LN*M*(75*DS4*M + 144*DS3*M*(-4*X + z^(-1)) + 4*DS2*(28*M^3 + 54*p2 + 216*M*X^2 + (45*M)/z^2 - (135*M*X)/z)) + (1080*S)/z^2 - (120*DS1*(-9 + M^2*z^2))/z^5) + (1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)))/z^6)*Vfun(-1/2*(z^4*(-2*DS1^3*LN*M + DS0*DS1*LN*M*(DS2 + DS1*(-3*X + z^(-1))) + DS0^2*LN*M*(DS3 + DS2*(-3*X + z^(-1))) - (2*DS0^3*(M - M*X*z + Phi*z^2))/z^3))/DS0^3))/(DS0^3*z^5);
end


    return co
end

function PhidotCoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 1/2;
 else 
    co = (-3*DS1^4*(597 - 4321*LN + 1710*LN^2)*M^2*z^6 + 9*DS0*DS1^2*M^2*z^5*(-3*DS2*(-159 + 887*LN + 380*LN^2)*z + 100*DS1*(3 - 16*LN + 4*(-3 + 19*LN)*X*z)) + 1800*DS0^4*(3 - 7*M^2*z^2 + 2*X*z*(6 + 5*M^2*z^2)) + DS0^2*M*z^4*(-6*DS2^2*(741 - 4963*LN + 855*LN^2)*M*z^2 - 15*DS1*M*z*(47*DS3*(3 - 19*LN)*z + 84*DS2*(3 - 16*LN + 4*(-3 + 19*LN)*X*z)) + 20*DS1^2*(-135*(-3 + 16*LN)*M*X*z + 54*(-3 + 19*LN)*p2*z^2 + 216*(-3 + 19*LN)*M*X^2*z^2 - 3*M*(45 + 23*M^2*z^2) + LN*M*(585 + 437*M^2*z^2))) + 5*DS0^3*(14040*S*z^4 - 240*DS1*z*(-18 + 5*M^2*z^2) + z^4*(3240*SZ*z + M*(-3*M*z*(25*DS4*(3 - 19*LN)*z + 48*DS3*(3 - 16*LN + 4*(-3 + 19*LN)*X*z)) + 4*DS2*(-135*(-3 + 16*LN)*M*X*z + 54*(-3 + 19*LN)*p2*z^2 + 216*(-3 + 19*LN)*M*X^2*z^2 - 3*M*(45 + 28*M^2*z^2) + LN*M*(585 + 532*M^2*z^2))))))/(2*(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(1080*S*z^4 - 120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z))))));
end


    return co
end

function PhidotCoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = z;
end


    return co
end

function PhidotCoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 0;
end


    return co
end

function PhidotSrc(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = (9*DS1^2*M - 9*DS0*DS2*M - 2*DS0^2*(2*M^3 - 9*p2 + 9*M*X^2))/(24*DS0^2);
 else 
    co = -1/8*((2*M*z^3*(4*DS1^3*(-1 + 3*LN)*z + DS0*DS1*(2*DS2*(1 - 3*LN)*z + 3*DS1*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)) + DS0^2*(2*DS3*(1 - 3*LN)*z + 3*DS2*(1 - 2*LN + 2*(-1 + 3*LN)*X*z))))/DS0^3 + (180*z*(2*DS1^3*(-1 + 4*LN)*M*z^3 + 2*DS0^3*(M - 2*M*X*z + 3*Phi*z^2 + PhiZ*z^3) + DS0*DS1*M*z^2*(DS2*(1 - 4*LN)*z + DS1*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)) + DS0^2*M*z^2*(DS3*(1 - 4*LN)*z + DS2*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)))*(-30*DS1^3*LN*M^2*z^5 + 30*DS0^3*(-3 - 6*X*z + M^2*z^2 - 3*X^2*z^2) - 9*DS0*DS1*z^2*(-6*DS2*LN*M^2*z^3 + 5*DS1*(2 - LN*M^2*z^2 + 2*LN*M^2*X*z^3)) + DS0^2*z*(-180*Sdot*z^3 + 20*DS1*(-9 - 9*X*z + 2*M^2*z^2) + 3*LN*M^2*z^3*(-4*DS3*z + 5*DS2*(-1 + 2*X*z)))))/(DS0^2*(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(1080*S*z^4 - 120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z)))))) + (3*M*z*(2*DS0^3 + 4*DS1^3*LN*z^3 + DS0*DS1*LN*z^2*(-2*DS2*z + DS1*(-3 + 6*X*z)) + DS0^2*LN*z^2*(-2*DS3*z + DS2*(-3 + 6*X*z)))*(-3*DS1^4*(199 - 1175*LN + 450*LN^2)*M^2*z^6 + 1800*DS0^4*(-3 - M^2*z^2 + 2*M^2*X*z^3) + 9*DS0*DS1^2*M^2*z^5*(-3*DS2*(-53 + 225*LN + 100*LN^2)*z + 100*DS1*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + DS0^2*M*z^4*(-6*DS2^2*(247 - 1325*LN + 225*LN^2)*M*z^2 - 15*DS1*M*z*(47*DS3*(1 - 5*LN)*z + 84*DS2*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + 20*DS1^2*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*z^2 + 216*(-1 + 5*LN)*M*X^2*z^2 + M*(-45 - 23*M^2*z^2 + 5*LN*(27 + 23*M^2*z^2)))) + 5*DS0^3*z^3*(-240*DS1*M^2 + 3240*S*z + z*(1080*SZ*z + M*(3*M*z*(25*DS4*(-1 + 5*LN)*z - 48*DS3*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + 4*DS2*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*z^2 + 216*(-1 + 5*LN)*M*X^2*z^2 + M*(-45 - 28*M^2*z^2 + 5*LN*(27 + 28*M^2*z^2))))))))/(DS0^3*(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(1080*S*z^4 - 120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z)))))) - 4*DV(-1/2*(z^4*(-2*DS1^3*LN*M + DS0*DS1*LN*M*(DS2 + DS1*(-3*X + z^(-1))) + DS0^2*LN*M*(DS3 + DS2*(-3*X + z^(-1))) - (2*DS0^3*(M - M*X*z + Phi*z^2))/z^3))/DS0^3))/z^3;
end


    return co
end

function ACoeff0(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 1;
 else 
    co = 6;
end


    return co
end

function ACoeff1(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = 6*z;
end


    return co
end

function ACoeff2(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = 0;
 else 
    co = z^2;
end


    return co
end

function ASrc(v::PDEVars{T}) where {T}

    @unpack Phi, S, Sdot, Phidot, A, PhiZ, SZ, SdotZ, PhidotZ, AZ, PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ, z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4 = v

if z == 0 
 co = -a4;
 else 
    co = (12 + (2*DS1*DS2*(7 - 12*LN)*M^2*z^5)/DS0^2 + (2*M^2*z^4*(DS3*(-7 + 12*LN)*z - 2*DS2*(5 - 6*LN + 2*(-7 + 12*LN)*X*z)))/DS0 + (3*z^2*(4*DS1^3*LN*M*z^3 + 2*DS0^3*(M - 2*Phidot*z^2) + DS0*DS1*LN*M*z^2*(-2*DS2*z + DS1*(-3 + 6*X*z)) + DS0^2*LN*M*z^2*(-2*DS3*z + DS2*(-3 + 6*X*z)))*(2*DS1^3*(-1 + 4*LN)*M*z^3 + 2*DS0^3*(M - 2*M*X*z + 3*Phi*z^2 + PhiZ*z^3) + DS0*DS1*M*z^2*(DS2*(1 - 4*LN)*z + DS1*(1 - 3*LN + 3*(-1 + 4*LN)*X*z)) + DS0^2*M*z^2*(DS3*(1 - 4*LN)*z + DS2*(1 - 3*LN + 3*(-1 + 4*LN)*X*z))))/DS0^6 - (2160*DS0*(-30*DS1^3*LN*M^2*z^5 + 30*DS0^3*(-3 - 6*X*z + M^2*z^2 - 3*X^2*z^2) - 9*DS0*DS1*z^2*(-6*DS2*LN*M^2*z^3 + 5*DS1*(2 - LN*M^2*z^2 + 2*LN*M^2*X*z^3)) + DS0^2*z*(-180*Sdot*z^3 + 20*DS1*(-9 - 9*X*z + 2*M^2*z^2) + 3*LN*M^2*z^3*(-4*DS3*z + 5*DS2*(-1 + 2*X*z))))*(-3*DS1^4*(199 - 1175*LN + 450*LN^2)*M^2*z^6 + 1800*DS0^4*(-3 - M^2*z^2 + 2*M^2*X*z^3) + 9*DS0*DS1^2*M^2*z^5*(-3*DS2*(-53 + 225*LN + 100*LN^2)*z + 100*DS1*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + DS0^2*M*z^4*(-6*DS2^2*(247 - 1325*LN + 225*LN^2)*M*z^2 - 15*DS1*M*z*(47*DS3*(1 - 5*LN)*z + 84*DS2*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + 20*DS1^2*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*z^2 + 216*(-1 + 5*LN)*M*X^2*z^2 + M*(-45 - 23*M^2*z^2 + 5*LN*(27 + 23*M^2*z^2)))) + 5*DS0^3*z^3*(-240*DS1*M^2 + 3240*S*z + z*(1080*SZ*z + M*(3*M*z*(25*DS4*(-1 + 5*LN)*z - 48*DS3*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)) + 4*DS2*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*z^2 + 216*(-1 + 5*LN)*M*X^2*z^2 + M*(-45 - 28*M^2*z^2 + 5*LN*(27 + 28*M^2*z^2))))))))/(3*DS1^4*(199 - 90*LN)*LN*M^2*z^6 - 9*DS0*DS1^2*LN*M^2*z^5*(3*DS2*(53 + 20*LN)*z - 100*DS1*(-1 + 4*X*z)) + 1800*DS0^4*(3 - M^2*z^2 + X*z*(3 + M^2*z^2)) + DS0^2*LN*M*z^4*(6*DS2^2*(247 - 45*LN)*M*z^2 + 20*DS1^2*(45*M - 135*M*X*z + 23*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) - 15*DS1*M*z*(-47*DS3*z + 84*DS2*(-1 + 4*X*z))) + 5*DS0^3*(1080*S*z^4 - 120*DS1*z*(-9 + M^2*z^2) + LN*M*z^4*(4*DS2*(45*M - 135*M*X*z + 28*M^3*z^2 + 54*p2*z^2 + 216*M*X^2*z^2) + 3*M*z*(25*DS4*z - 48*DS3*(-1 + 4*X*z)))))^2 - 8*Vfun(-1/2*(z^4*(-2*DS1^3*LN*M + DS0*DS1*LN*M*(DS2 + DS1*(-3*X + z^(-1))) + DS0^2*LN*M*(DS3 + DS2*(-3*X + z^(-1))) - (2*DS0^3*(M - M*X*z + Phi*z^2))/z^3))/DS0^3))/(6*z^4);
end


    return co
end


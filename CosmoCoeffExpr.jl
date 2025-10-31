function DS0(t)

res = 0.857763885*Power(E,t/2.)*Sqrt(Cosh(1 - t)); 

return res;
end

function DS1(t)

res = (0.157776849*Power(E,(3*t)/2.))/Sqrt(Cosh(1 - t)); 

return res;
end

function DS2(t)

res = (0.02902143*Power(E,t/2.)*(14.7781122 + 1. *Power(E,2*t)))/Power(Cosh(1 - t),1.5); 

return res;
end

function DS3(t)

res = (1.165822 + 0.07888842*Power(E,2*t) + 0.005338194*Power(E,4*t))/(Power(E,t/2.)*Power(Cosh(1 - t),2.5)); 

return res;
end

function DS4(t)

res = (3.169033 - 0.214441*Power(E,2*t) + 0.02902143*Power(E,4*t) + 0.000981906*Power(E,6*t))/(Power(E,(3*t)/2.)*Power(Cosh(1 - t),3.5)); 

return res;
end

function DS5(t)

res = (8.6 - 4.08*Power(E,2*t) + 0.355*Power(E,4*t) + 0.00534*Power(E,6*t) + 0.000181*Power(E,8*t))/(Power(E,(5*t)/2.)*Power(Cosh(1 - t),4.5)); 

return res;
end

function DS6(t)

res = (23.416 - 33.27*Power(E,2*t) + 6.648*Power(E,4*t) - 0.13785*Power(E,6*t) + 0.0014729*Power(E,8*t) + 0.00003322*Power(E,10*t))/(Power(E,(7*t)/2.)*Power(Cosh(1 - t),5.5)); 

return res;
end

function SCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (72*Power(DS0(t),6) + Power(2*z*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3))*Power(DS0(t),3) + 2*(-1 + 4*LN)*M*Power(z,4)*Power(DS1(t),3) + M*Power(z,3)*DS0(t)*DS1(t)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS1(t) + (1 - 4*LN)*z*DS2(t)) + M*Power(z,3)*Power(DS0(t),2)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS2(t) + (1 - 4*LN)*z*DS3(t)),2))/(6. *Power(z,2)*Power(DS0(t),6)); 

return co;
end

function SCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 8/z; 

return co;
end

function SSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (Power(2*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3))*Power(DS0(t),3) + 2*(-1 + 4*LN)*M*Power(z,3)*Power(DS1(t),3) + M*Power(z,2)*DS0(t)*DS1(t)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS1(t) + (1 - 4*LN)*z*DS2(t)) + M*Power(z,2)*Power(DS0(t),2)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS2(t) + (1 - 4*LN)*z*DS3(t)),2)*(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(-120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t))))) + 6*M*Power(DS0(t),6)*(3600*M*(-1 + 3*X*z)*Power(DS0(t),4) - 3*(2369 - 7950*LN + 2700*Power(LN,2))*M*Power(z,4)*Power(DS1(t),4) + 9*M*Power(z,3)*DS0(t)*Power(DS1(t),2)*(100*(9 - 20*LN + 4*(-11 + 30*LN)*X*z)*DS1(t) - 3*(-543 + 1150*LN + 600*Power(LN,2))*z*DS2(t)) + Power(z,2)*Power(DS0(t),2)*(20*(-135*(-9 + 20*LN)*M*X*z + 54*(-11 + 30*LN)*p2*Power(z,2) + 216*(-11 + 30*LN)*M*Power(X,2)*Power(z,2) + M*(-315 - 253*Power(M,2)*Power(z,2) + 30*LN*(18 + 23*Power(M,2)*Power(z,2))))*Power(DS1(t),2) - 6*(2807 - 8400*LN + 1350*Power(LN,2))*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(9 - 20*LN + 4*(-11 + 30*LN)*X*z)*DS2(t) + 47*(11 - 30*LN)*z*DS3(t))) + 5*z*Power(DS0(t),3)*(-720*M*DS1(t) + z*(4*(-135*(-9 + 20*LN)*M*X*z + 54*(-11 + 30*LN)*p2*Power(z,2) + 216*(-11 + 30*LN)*M*Power(X,2)*Power(z,2) + 60*LN*M*(9 + 14*Power(M,2)*Power(z,2)) - 7*M*(45 + 44*Power(M,2)*Power(z,2)))*DS2(t) + 3*M*z*(-48*(9 - 20*LN + 4*(-11 + 30*LN)*X*z)*DS3(t) + 25*(-11 + 30*LN)*z*DS4(t))))))/(32400. *Power(z,4)*Power(DS0(t),9)); 

return co;
end

function SdotCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (2*(1800*(-2*Power(M,2)*z + 3*X*(1 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) - 3*(199 - 1374*LN + 540*Power(LN,2))*Power(M,2)*Power(z,5)*Power(DS1(t),4) + 9*Power(M,2)*Power(z,4)*DS0(t)*Power(DS1(t),2)*(100*(1 - 5*LN + 4*(-1 + 6*LN)*X*z)*DS1(t) - 3*(-53 + 278*LN + 120*Power(LN,2))*z*DS2(t)) + M*Power(z,3)*Power(DS0(t),2)*(20*(-135*(-1 + 5*LN)*M*X*z + 54*(-1 + 6*LN)*p2*Power(z,2) + 216*(-1 + 6*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 23*Power(M,2)*Power(z,2) + 6*LN*(30 + 23*Power(M,2)*Power(z,2))))*Power(DS1(t),2) - 6*(247 - 1572*LN + 270*Power(LN,2))*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(1 - 5*LN + 4*(-1 + 6*LN)*X*z)*DS2(t) + 47*(1 - 6*LN)*z*DS3(t))) + 5*Power(DS0(t),3)*(4320*S*Power(z,3) - 360*(-3 + Power(M,2)*Power(z,2))*DS1(t) + Power(z,3)*(1080*SZ*z + M*(4*(-135*(-1 + 5*LN)*M*X*z + 54*(-1 + 6*LN)*p2*Power(z,2) + 216*(-1 + 6*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 28*Power(M,2)*Power(z,2) + 12*LN*(15 + 14*Power(M,2)*Power(z,2))))*DS2(t) + 3*M*z*(-48*(1 - 5*LN + 4*(-1 + 6*LN)*X*z)*DS3(t) + 25*(-1 + 6*LN)*z*DS4(t)))))))/(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(1080*S*Power(z,4) - 120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t))))); 

return co;
end

function SdotCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 1; 

return co;
end

function SdotSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = -0.0001234567901234568*(45*DS0(t)*(180*(1 + X*z)*Power(DS0(t),3) + 30*(1 - 3*LN)*Power(M,2)*Power(z,5)*Power(DS1(t),3) - 9*Power(M,2)*Power(z,4)*DS0(t)*DS1(t)*(5*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)*DS1(t) + 6*(1 - 3*LN)*z*DS2(t)) + Power(DS0(t),2)*(20*z*(9 + 2*Power(M,2)*Power(z,2))*DS1(t) + 3*Power(M,2)*Power(z,4)*(5*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)*DS2(t) + 4*(1 - 3*LN)*z*DS3(t)))) + (90*DS0(t)*(30*(-3 - 6*X*z + Power(M,2)*Power(z,2) - 3*Power(X,2)*Power(z,2))*Power(DS0(t),3) - 30*LN*Power(M,2)*Power(z,5)*Power(DS1(t),3) - 9*Power(z,2)*DS0(t)*DS1(t)*(5*(2 - LN*Power(M,2)*Power(z,2) + 2*LN*Power(M,2)*X*Power(z,3))*DS1(t) - 6*LN*Power(M,2)*Power(z,3)*DS2(t)) + Power(DS0(t),2)*(20*z*(-9 - 9*X*z + 2*Power(M,2)*Power(z,2))*DS1(t) + 3*LN*Power(M,2)*Power(z,4)*(5*(-1 + 2*X*z)*DS2(t) - 4*z*DS3(t))))*(1800*(-3 - Power(M,2)*Power(z,2) + 2*Power(M,2)*X*Power(z,3))*Power(DS0(t),4) - 3*(199 - 1175*LN + 450*Power(LN,2))*Power(M,2)*Power(z,6)*Power(DS1(t),4) + 9*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(100*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS1(t) - 3*(-53 + 225*LN + 100*Power(LN,2))*z*DS2(t)) + M*Power(z,4)*Power(DS0(t),2)*(20*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*Power(z,2) + 216*(-1 + 5*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 23*Power(M,2)*Power(z,2) + 5*LN*(27 + 23*Power(M,2)*Power(z,2))))*Power(DS1(t),2) - 6*(247 - 1325*LN + 225*Power(LN,2))*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS2(t) + 47*(1 - 5*LN)*z*DS3(t))) + 5*Power(z,3)*Power(DS0(t),3)*(3240*S*z - 240*Power(M,2)*DS1(t) + z*(1080*SZ*z + M*(4*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*Power(z,2) + 216*(-1 + 5*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 28*Power(M,2)*Power(z,2) + 5*LN*(27 + 28*Power(M,2)*Power(z,2))))*DS2(t) + 3*M*z*(-48*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS3(t) + 25*(-1 + 5*LN)*z*DS4(t)))))))/(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(1080*S*Power(z,4) - 120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t))))) + Power(z,6)*((1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4))/Power(z,6) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(DS1(t),4) - 9*LN*Power(M,2)*DS0(t)*Power(DS1(t),2)*(100*(-4*X + 1/z)*DS1(t) + 3*(53 + 20*LN)*DS2(t)) + LN*M*Power(DS0(t),2)*(20*(23*Power(M,3) + 54*p2 + 216*M*Power(X,2) + (45*M)/Power(z,2) - (135*M*X)/z)*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(DS2(t),2) + 15*M*DS1(t)*(84*(-4*X + 1/z)*DS2(t) + 47*DS3(t))) + 5*Power(DS0(t),3)*((1080*S)/Power(z,2) - (120*(-9 + Power(M,2)*Power(z,2))*DS1(t))/Power(z,5) + LN*M*(4*(28*Power(M,3) + 54*p2 + 216*M*Power(X,2) + (45*M)/Power(z,2) - (135*M*X)/z)*DS2(t) + 144*M*(-4*X + 1/z)*DS3(t) + 75*M*DS4(t))))*Vfun(-0.5*(Power(z,4)*((-2*(M - M*X*z + Phi*Power(z,2))*Power(DS0(t),3))/Power(z,3) - 2*LN*M*Power(DS1(t),3) + LN*M*DS0(t)*DS1(t)*((-3*X + 1/z)*DS1(t) + DS2(t)) + LN*M*Power(DS0(t),2)*((-3*X + 1/z)*DS2(t) + DS3(t))))/Power(DS0(t),3)))/(Power(z,5)*Power(DS0(t),3)); 

return co;
end

function PhidotCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (1800*(3 - 7*Power(M,2)*Power(z,2) + 2*X*z*(6 + 5*Power(M,2)*Power(z,2)))*Power(DS0(t),4) - 3*(597 - 4321*LN + 1710*Power(LN,2))*Power(M,2)*Power(z,6)*Power(DS1(t),4) + 9*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(100*(3 - 16*LN + 4*(-3 + 19*LN)*X*z)*DS1(t) - 3*(-159 + 887*LN + 380*Power(LN,2))*z*DS2(t)) + M*Power(z,4)*Power(DS0(t),2)*(20*(-135*(-3 + 16*LN)*M*X*z + 54*(-3 + 19*LN)*p2*Power(z,2) + 216*(-3 + 19*LN)*M*Power(X,2)*Power(z,2) - 3*M*(45 + 23*Power(M,2)*Power(z,2)) + LN*M*(585 + 437*Power(M,2)*Power(z,2)))*Power(DS1(t),2) - 6*(741 - 4963*LN + 855*Power(LN,2))*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(3 - 16*LN + 4*(-3 + 19*LN)*X*z)*DS2(t) + 47*(3 - 19*LN)*z*DS3(t))) + 5*Power(DS0(t),3)*(14040*S*Power(z,4) - 240*z*(-18 + 5*Power(M,2)*Power(z,2))*DS1(t) + Power(z,4)*(3240*SZ*z + M*(4*(-135*(-3 + 16*LN)*M*X*z + 54*(-3 + 19*LN)*p2*Power(z,2) + 216*(-3 + 19*LN)*M*Power(X,2)*Power(z,2) - 3*M*(45 + 28*Power(M,2)*Power(z,2)) + LN*M*(585 + 532*Power(M,2)*Power(z,2)))*DS2(t) - 3*M*z*(48*(3 - 16*LN + 4*(-3 + 19*LN)*X*z)*DS3(t) + 25*(3 - 19*LN)*z*DS4(t))))))/(2. *z*(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(1080*S*Power(z,4) - 120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t)))))); 

return co;
end

function PhidotCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 1; 

return co;
end

function PhidotSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = -0.125*((2*M*Power(z,3)*(4*(-1 + 3*LN)*z*Power(DS1(t),3) + DS0(t)*DS1(t)*(3*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)*DS1(t) + 2*(1 - 3*LN)*z*DS2(t)) + Power(DS0(t),2)*(3*(1 - 2*LN + 2*(-1 + 3*LN)*X*z)*DS2(t) + 2*(1 - 3*LN)*z*DS3(t))))/Power(DS0(t),3) + (180*z*(2*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3))*Power(DS0(t),3) + 2*(-1 + 4*LN)*M*Power(z,3)*Power(DS1(t),3) + M*Power(z,2)*DS0(t)*DS1(t)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS1(t) + (1 - 4*LN)*z*DS2(t)) + M*Power(z,2)*Power(DS0(t),2)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS2(t) + (1 - 4*LN)*z*DS3(t)))*(30*(-3 - 6*X*z + Power(M,2)*Power(z,2) - 3*Power(X,2)*Power(z,2))*Power(DS0(t),3) - 30*LN*Power(M,2)*Power(z,5)*Power(DS1(t),3) - 9*Power(z,2)*DS0(t)*DS1(t)*(5*(2 - LN*Power(M,2)*Power(z,2) + 2*LN*Power(M,2)*X*Power(z,3))*DS1(t) - 6*LN*Power(M,2)*Power(z,3)*DS2(t)) + z*Power(DS0(t),2)*(-180*Sdot*Power(z,3) + 20*(-9 - 9*X*z + 2*Power(M,2)*Power(z,2))*DS1(t) + 3*LN*Power(M,2)*Power(z,3)*(5*(-1 + 2*X*z)*DS2(t) - 4*z*DS3(t)))))/(Power(DS0(t),2)*(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(1080*S*Power(z,4) - 120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t)))))) + (3*M*z*(2*Power(DS0(t),3) + 4*LN*Power(z,3)*Power(DS1(t),3) + LN*Power(z,2)*DS0(t)*DS1(t)*((-3 + 6*X*z)*DS1(t) - 2*z*DS2(t)) + LN*Power(z,2)*Power(DS0(t),2)*((-3 + 6*X*z)*DS2(t) - 2*z*DS3(t)))*(1800*(-3 - Power(M,2)*Power(z,2) + 2*Power(M,2)*X*Power(z,3))*Power(DS0(t),4) - 3*(199 - 1175*LN + 450*Power(LN,2))*Power(M,2)*Power(z,6)*Power(DS1(t),4) + 9*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(100*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS1(t) - 3*(-53 + 225*LN + 100*Power(LN,2))*z*DS2(t)) + M*Power(z,4)*Power(DS0(t),2)*(20*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*Power(z,2) + 216*(-1 + 5*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 23*Power(M,2)*Power(z,2) + 5*LN*(27 + 23*Power(M,2)*Power(z,2))))*Power(DS1(t),2) - 6*(247 - 1325*LN + 225*Power(LN,2))*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS2(t) + 47*(1 - 5*LN)*z*DS3(t))) + 5*Power(z,3)*Power(DS0(t),3)*(3240*S*z - 240*Power(M,2)*DS1(t) + z*(1080*SZ*z + M*(4*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*Power(z,2) + 216*(-1 + 5*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 28*Power(M,2)*Power(z,2) + 5*LN*(27 + 28*Power(M,2)*Power(z,2))))*DS2(t) + 3*M*z*(-48*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS3(t) + 25*(-1 + 5*LN)*z*DS4(t)))))))/(Power(DS0(t),3)*(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(1080*S*Power(z,4) - 120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t)))))) - 4*DV(-0.5*(Power(z,4)*((-2*(M - M*X*z + Phi*Power(z,2))*Power(DS0(t),3))/Power(z,3) - 2*LN*M*Power(DS1(t),3) + LN*M*DS0(t)*DS1(t)*((-3*X + 1/z)*DS1(t) + DS2(t)) + LN*M*Power(DS0(t),2)*((-3*X + 1/z)*DS2(t) + DS3(t))))/Power(DS0(t),3)))/Power(z,4); 

return co;
end

function ACoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 6/Power(z,2); 

return co;
end

function ACoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 6/z; 

return co;
end

function ASrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (12 + (2*(7 - 12*LN)*Power(M,2)*Power(z,5)*DS1(t)*DS2(t))/Power(DS0(t),2) + (2*Power(M,2)*Power(z,4)*(-2*(5 - 6*LN + 2*(-7 + 12*LN)*X*z)*DS2(t) + (-7 + 12*LN)*z*DS3(t)))/DS0(t) + (3*Power(z,2)*(2*(M - 2*Phidot*Power(z,2))*Power(DS0(t),3) + 4*LN*M*Power(z,3)*Power(DS1(t),3) + LN*M*Power(z,2)*DS0(t)*DS1(t)*((-3 + 6*X*z)*DS1(t) - 2*z*DS2(t)) + LN*M*Power(z,2)*Power(DS0(t),2)*((-3 + 6*X*z)*DS2(t) - 2*z*DS3(t)))*(2*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3))*Power(DS0(t),3) + 2*(-1 + 4*LN)*M*Power(z,3)*Power(DS1(t),3) + M*Power(z,2)*DS0(t)*DS1(t)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS1(t) + (1 - 4*LN)*z*DS2(t)) + M*Power(z,2)*Power(DS0(t),2)*((1 - 3*LN + 3*(-1 + 4*LN)*X*z)*DS2(t) + (1 - 4*LN)*z*DS3(t))))/Power(DS0(t),6) - (2160*DS0(t)*(30*(-3 - 6*X*z + Power(M,2)*Power(z,2) - 3*Power(X,2)*Power(z,2))*Power(DS0(t),3) - 30*LN*Power(M,2)*Power(z,5)*Power(DS1(t),3) - 9*Power(z,2)*DS0(t)*DS1(t)*(5*(2 - LN*Power(M,2)*Power(z,2) + 2*LN*Power(M,2)*X*Power(z,3))*DS1(t) - 6*LN*Power(M,2)*Power(z,3)*DS2(t)) + z*Power(DS0(t),2)*(-180*Sdot*Power(z,3) + 20*(-9 - 9*X*z + 2*Power(M,2)*Power(z,2))*DS1(t) + 3*LN*Power(M,2)*Power(z,3)*(5*(-1 + 2*X*z)*DS2(t) - 4*z*DS3(t))))*(1800*(-3 - Power(M,2)*Power(z,2) + 2*Power(M,2)*X*Power(z,3))*Power(DS0(t),4) - 3*(199 - 1175*LN + 450*Power(LN,2))*Power(M,2)*Power(z,6)*Power(DS1(t),4) + 9*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(100*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS1(t) - 3*(-53 + 225*LN + 100*Power(LN,2))*z*DS2(t)) + M*Power(z,4)*Power(DS0(t),2)*(20*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*Power(z,2) + 216*(-1 + 5*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 23*Power(M,2)*Power(z,2) + 5*LN*(27 + 23*Power(M,2)*Power(z,2))))*Power(DS1(t),2) - 6*(247 - 1325*LN + 225*Power(LN,2))*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS2(t) + 47*(1 - 5*LN)*z*DS3(t))) + 5*Power(z,3)*Power(DS0(t),3)*(3240*S*z - 240*Power(M,2)*DS1(t) + z*(1080*SZ*z + M*(4*(-135*(-1 + 4*LN)*M*X*z + 54*(-1 + 5*LN)*p2*Power(z,2) + 216*(-1 + 5*LN)*M*Power(X,2)*Power(z,2) + M*(-45 - 28*Power(M,2)*Power(z,2) + 5*LN*(27 + 28*Power(M,2)*Power(z,2))))*DS2(t) + 3*M*z*(-48*(1 - 4*LN + 4*(-1 + 5*LN)*X*z)*DS3(t) + 25*(-1 + 5*LN)*z*DS4(t)))))))/Power(1800*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Power(DS0(t),4) + 3*(199 - 90*LN)*LN*Power(M,2)*Power(z,6)*Power(DS1(t),4) - 9*LN*Power(M,2)*Power(z,5)*DS0(t)*Power(DS1(t),2)*(-100*(-1 + 4*X*z)*DS1(t) + 3*(53 + 20*LN)*z*DS2(t)) + LN*M*Power(z,4)*Power(DS0(t),2)*(20*(45*M - 135*M*X*z + 23*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*Power(DS1(t),2) + 6*(247 - 45*LN)*M*Power(z,2)*Power(DS2(t),2) - 15*M*z*DS1(t)*(84*(-1 + 4*X*z)*DS2(t) - 47*z*DS3(t))) + 5*Power(DS0(t),3)*(1080*S*Power(z,4) - 120*z*(-9 + Power(M,2)*Power(z,2))*DS1(t) + LN*M*Power(z,4)*(4*(45*M - 135*M*X*z + 28*Power(M,3)*Power(z,2) + 54*p2*Power(z,2) + 216*M*Power(X,2)*Power(z,2))*DS2(t) + 3*M*z*(-48*(-1 + 4*X*z)*DS3(t) + 25*z*DS4(t)))),2) - 8*Vfun(-0.5*(Power(z,4)*((-2*(M - M*X*z + Phi*Power(z,2))*Power(DS0(t),3))/Power(z,3) - 2*LN*M*Power(DS1(t),3) + LN*M*DS0(t)*DS1(t)*((-3*X + 1/z)*DS1(t) + DS2(t)) + LN*M*Power(DS0(t),2)*((-3*X + 1/z)*DS2(t) + DS3(t))))/Power(DS0(t),3)))/(6. *Power(z,6)); 

return co;
end


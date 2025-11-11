function SCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (2*(18 + Power(M,2)*Power(z,2) + 4*Power(M,2)*Power(X,2)*Power(z,4) + 2*M*PhiZ*Power(z,5) + 9*Power(Phi,2)*Power(z,6) + Power(PhiZ,2)*Power(z,8) - 4*M*X*Power(z,3)*(M + PhiZ*Power(z,3)) + 6*Phi*Power(z,4)*(M - 2*M*X*z + PhiZ*Power(z,3))))/3.; 

return co;
end

function SCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = 8*z; 

return co;
end

function SCoeff2(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = Power(z,2); 

return co;
end

function SSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (2*(-Power(M,4) + 6*M*PhiZ*z - 2*Power(M,3)*PhiZ*Power(z,3) + 3*Power(PhiZ,2)*Power(z,4) - Power(M,2)*Power(PhiZ,2)*Power(z,6) + 4*Power(M,2)*Power(X,3)*z*(3 + Power(M,2)*Power(z,2)) - 4*M*Power(X,2)*Power(z,2)*(2*Power(M,3) + PhiZ*z*(3 + Power(M,2)*Power(z,2))) + 9*Power(Phi,2)*Power(z,2)*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2))) - 6*Phi*(-M + 2*M*X*z - PhiZ*Power(z,3))*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2))) + X*z*(5*Power(M,4) + 6*M*PhiZ*z*(-1 + Power(M,2)*Power(z,2)) + Power(PhiZ,2)*Power(z,4)*(3 + Power(M,2)*Power(z,2)))))/9.; 

return co;
end

function SdotCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (2*Power(z,5)*(-2*Power(M,2)*z + 12*S*Power(z,3) + 3*SZ*Power(z,4) + 3*X*(1 + Power(M,2)*Power(z,2))))/(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2))); 

return co;
end

function SdotCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = Power(z,5); 

return co;
end

function SdotCoeff2(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = 0; 

return co;
end

function SdotSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (18*Power(M,2)*Power(X,3)*Power(z,5) + 3*(-3 + Power(M,2)*Power(z,2))*(6 + Power(M,2)*Power(z,2) - 3*SZ*Power(z,5)) + 9*Power(X,2)*Power(z,2)*(-6 + 2*Power(M,2)*Power(z,2) + 3*SZ*Power(z,5)) + 6*X*z*(-18 - Power(M,4)*Power(z,4) + 9*SZ*Power(z,5)) - 18*Power(S,2)*Power(z,8)*Vfun(z*(M - M*X*z + Phi*Power(z,2))) - 2*Power(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)),2)*Vfun(z*(M - M*X*z + Phi*Power(z,2))) - 3*S*Power(z,4)*(9*(-2 - 5*X*z + Power(M,2)*Power(z,2) - 3*Power(X,2)*Power(z,2)) + 4*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))*Vfun(z*(M - M*X*z + Phi*Power(z,2)))))/(9. *(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))); 

return co;
end

function PhidotCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (3 - 7*Power(M,2)*Power(z,2) + 39*S*Power(z,4) + 9*SZ*Power(z,5) + 2*X*z*(6 + 5*Power(M,2)*Power(z,2)))/(6 - 2*Power(M,2)*Power(z,2) + 6*S*Power(z,4) + 2*X*z*(3 + Power(M,2)*Power(z,2))); 

return co;
end

function PhidotCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = z; 

return co;
end

function PhidotCoeff2(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = 0; 

return co;
end

function PhidotSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (18*M*z - 27*M*Power(X,2)*Power(z,3) + 9*PhiZ*Power(z,4) - 18*M*Power(X,3)*Power(z,4) - 27*M*S*Power(z,5) + 18*PhiZ*X*Power(z,5) - 3*Power(M,2)*PhiZ*Power(z,6) - 9*M*SZ*Power(z,6) + 9*PhiZ*Power(X,2)*Power(z,6) + 18*Sdot*Power(z,5)*(M - 2*M*X*z + PhiZ*Power(z,3)) + 9*Phi*Power(z,3)*(3 + 6*X*z - Power(M,2)*Power(z,2) + 3*Power(X,2)*Power(z,2) + 6*Sdot*Power(z,4)) + 6*DV(z*(M - M*X*z + Phi*Power(z,2))) + 6*X*z*DV(z*(M - M*X*z + Phi*Power(z,2))) - 2*Power(M,2)*Power(z,2)*DV(z*(M - M*X*z + Phi*Power(z,2))) + 2*Power(M,2)*X*Power(z,3)*DV(z*(M - M*X*z + Phi*Power(z,2))) + 6*S*Power(z,4)*DV(z*(M - M*X*z + Phi*Power(z,2))))/(4. *Power(z,3)*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))); 

return co;
end

function ACoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = 6; 

return co;
end

function ACoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = 6*z; 

return co;
end

function ACoeff2(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = Power(z,2); 

return co;
end

function ASrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

co = (2*(3 + 3*Power(z,2)*(M - 2*Phidot*Power(z,2))*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3)) + (9*(3 + 6*X*z - Power(M,2)*Power(z,2) + 3*Power(X,2)*Power(z,2) + 6*Sdot*Power(z,4))*(-3 - Power(M,2)*Power(z,2) + 2*Power(M,2)*X*Power(z,3) + 9*S*Power(z,4) + 3*SZ*Power(z,5)))/Power(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)),2) - 2*Vfun(z*(M - M*X*z + Phi*Power(z,2)))))/(3. *Power(z,4)); 

return co;
end


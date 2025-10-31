function SCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (72 + 4*Power(z,2)*Power(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3),2))/(6. *Power(z,2)); 

return co;
end

function SCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 8/z; 

return co;
end

function SSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (2*(3*Power(M,2)*(-1 + 3*X*z) + Power(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3),2)*(3 - Power(M,2)*Power(z,2) + X*z*(3 + Power(M,2)*Power(z,2)))))/(9. *Power(z,4)); 

return co;
end

function SdotCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (-4*Power(M,2)*z + 24*S*Power(z,3) + 6*SZ*Power(z,4) + 6*X*(1 + Power(M,2)*Power(z,2)))/(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2))); 

return co;
end

function SdotCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 1; 

return co;
end

function SdotSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = -0.1111111111111111*(9 + 9*X*z + (3*(-3 - 6*X*z + Power(M,2)*Power(z,2) - 3*Power(X,2)*Power(z,2))*(-3 - Power(M,2)*Power(z,2) + 2*Power(M,2)*X*Power(z,3) + 9*S*Power(z,4) + 3*SZ*Power(z,5)))/(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2))) + 2*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))*Vfun(z*(M - M*X*z + Phi*Power(z,2))))/Power(z,5); 

return co;
end

function PhidotCoeff0(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (3 - 7*Power(M,2)*Power(z,2) + 39*S*Power(z,4) + 9*SZ*Power(z,5) + 2*X*z*(6 + 5*Power(M,2)*Power(z,2)))/(2. *z*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))); 

return co;
end

function PhidotCoeff1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = 1; 

return co;
end

function PhidotSrc(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4)

co = (3*z*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3))*(3 + 6*X*z - Power(M,2)*Power(z,2) + 3*Power(X,2)*Power(z,2) + 6*Sdot*Power(z,4)) + 3*M*z*(3 + Power(M,2)*Power(z,2) - 2*Power(M,2)*X*Power(z,3) - 9*S*Power(z,4) - 3*SZ*Power(z,5)) + 2*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))*DV(z*(M - M*X*z + Phi*Power(z,2))))/(4. *Power(z,4)*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))); 

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

co = (2*(3 + 3*Power(z,2)*(M - 2*Phidot*Power(z,2))*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3)) + (9*(3 + 6*X*z - Power(M,2)*Power(z,2) + 3*Power(X,2)*Power(z,2) + 6*Sdot*Power(z,4))*(-3 - Power(M,2)*Power(z,2) + 2*Power(M,2)*X*Power(z,3) + 9*S*Power(z,4) + 3*SZ*Power(z,5)))/Power(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)),2) - 2*Vfun(z*(M - M*X*z + Phi*Power(z,2)))))/(3. *Power(z,6)); 

return co;
end


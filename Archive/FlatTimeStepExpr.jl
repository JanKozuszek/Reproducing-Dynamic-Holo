function DtX(S,SZ, Sdot,SdotZ, Phidot, A, AZ, z, LN, X, p2, t)

XPrime = (Power(z,2)*(2*Power(M,4) + 72*Sdot - 9*AZ*z + 18*SdotZ*z - 2*Power(M,4)*X*z + 108*Sdot*X*z - 6*Power(M,2)*S*Power(z,2) - 24*Power(M,2)*Sdot*Power(z,2) - 18*AZ*X*Power(z,2) + 36*SdotZ*X*Power(z,2) + 36*Sdot*Power(X,2)*Power(z,2) + 3*AZ*Power(M,2)*Power(z,3) - 12*Power(M,2)*SdotZ*Power(z,3) - 9*AZ*Power(X,2)*Power(z,3) + 18*SdotZ*Power(X,2)*Power(z,3) - 18*AZ*Sdot*Power(z,5) + 6*A*(-6 - 9*X*z + Power(M,2)*Power(z,2) - 3*Power(X,2)*Power(z,2) + 3*SdotZ*Power(z,5)) + 8*M*Phidot*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2))) - 8*Power(Phidot,2)*Power(z,2)*(3 - Power(M,2)*Power(z,2) + 3*S*Power(z,4) + X*z*(3 + Power(M,2)*Power(z,2)))))/(36. *(-1 - X*z + 2*Sdot*Power(z,4) + SdotZ*Power(z,5)));

return XPrime;
end;

function DtPhi(Phi, PhiZ, Phidot, A, z, LN, X, XPrime, t)

timeder = (-2*Power(M,3) + 9*Phi + 6*Phidot - 9*M*Power(X,2) + 3*PhiZ*z + 4*Power(M,3)*X*z + 18*Phi*X*z - 6*M*Power(X,3)*z + 12*M*X*XPrime*z - 6*Power(M,2)*Phi*Power(z,2) + 6*PhiZ*X*Power(z,2) + 9*Phi*Power(X,2)*Power(z,2) - 18*Phi*XPrime*Power(z,2) - 2*Power(M,2)*PhiZ*Power(z,3) + 3*PhiZ*Power(X,2)*Power(z,3) - 6*PhiZ*XPrime*Power(z,3) + 3*A*Power(z,2)*(M - 2*M*X*z + 3*Phi*Power(z,2) + PhiZ*Power(z,3)))/(6. *z);

return timeder;
end;

function Dta4(X, a4, p2, XPrime, p2Prime, t)

a4Prime = (2*(-18*M*p2Prime + 36*Power(M,2)*X*XPrime))/27.;

return a4Prime;
end;


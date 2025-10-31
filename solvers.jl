# Functions for computing equations etc. Order of variables: Phi, S, Sdot, Phidot, A.
# X = ξ. XPrime = ξ'. In general, spatial derivatives are denoted by adding a Z at the end, i.e. AZZ = ∂_zz(A).
# We need the top couple of functions due to how Mathematica's CForm treats functions, e.g. x^y gets exported as Power(x,y).
# The data needed at the start of a timestep: a profile for Φ (which includes ϕ_2), values for a_4 and ξ.

import LinearAlgebra as linalg
import Plots as plt
import Polynomials as poly
import NonlinearSolve as nonlin

using Format
using .Threads
using ProgressMeter


function Power(x,y)
    return(x^y)
end

function Cosh(x)
    return cosh(x)
end

function Sech(x) 
    return sech(x)
end

function Sqrt(x)
    return sqrt(x)
end

function Unroll(arr)
    #Takes an array subdivided into Ndom domains, returns a flat vector of length N = Ndom*Npts.
    unrolled = zeros(T,N);
    for ii in 1:Ndom
        unrolled[(ii-1)*Npts+1:ii*Npts] = arr[ii,:];
    end
    return unrolled
end;

function ComputeDerivatives(Var)
    #Takes a full state vector of dimensions NVar × N, returns all first and second z-derivatives.
    VarZ  = zeros(T,NVar,N);
    VarZZ = zeros(T,NVar,N);

    for ii in 1:NVar
        Threads.@threads for jj in 1:Ndom
            locrange = (jj - 1)*Npts+1:jj*Npts;
            varhere = Var[ii, locrange];
            VarZ[ii, locrange]  = DiffMats[jj,:,:]*varhere;
            VarZZ[ii, locrange] = DiffMats2[jj,:,:]*varhere;
        end
    end

    return([VarZ,VarZZ])
end

function ComputeSingleDerivative(Vec)
    #Takes a single variable vector, returns its first derivative only.
    VecZ = zeros(T,N)

    Threads.@threads for jj in 1:Ndom
        locrange = (jj - 1)*Npts+1:jj*Npts;
        varhere = Vec[locrange];
        VecZ[locrange] = DiffMats[jj,:,:]*varhere;
    end
    return VecZ
    
end


function LinearSolveODE(Var, EqNum, a4, X, t)
    # Solve an ODE by direct linear inversion. Grid is subdivided into Ndom domains of Npts each.
    # This is a very complicated function. Maybe better split into bits.
    degreelist = [2,1,1,2];
    regBC = [true, false, true, false];
    varnamelist = ["S","Sdot","Phidot","A"];
    degree = degreelist[EqNum];
    VarZ, VarZZ = ComputeDerivatives(Var);
    p2 = Var[1,1];
    sol = zeros(T,N);
    solParticular = zeros(T,Ndom,Npts);
    solHomogeneous = zeros(T,Ndom,Npts);
    srcPrescribed = T(0.);

    if degree == 2
        solHomogeneous2 = zeros(T,Ndom,Npts);

        co1_functionName = string(varnamelist[EqNum],"Coeff1");
        co1_function = getfield(Main, Symbol(co1_functionName));
    end

    co0_functionName = string(varnamelist[EqNum],"Coeff0");
    co0_function = getfield(Main, Symbol(co0_functionName));

    src_functionName = string(varnamelist[EqNum],"Src");
    src_function = getfield(Main, Symbol(src_functionName));

    #Basic parallelization here:
    Threads.@threads for domind in 1:Ndom

        if degree == 2
            F1Mat = zeros(T,Npts,Npts);
        end
        F0Mat = zeros(T,Npts,Npts);
        SRCVec = zeros(T,Npts);
        HomoSRCVec = zeros(T, Npts);
        HomoSRCVec[1] = 1;

        for ptind in 1:Npts
            fullind = (domind - 1)*Npts + ptind;
            Phi, S, Sdot, Phidot, A = Var[1:NVar,fullind];
            PhiZ, SZ, SdotZ, PhidotZ, AZ = VarZ[1:NVar,fullind]
            PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ = VarZZ[1:NVar,fullind];
            z = grid[fullind]; LN = -log(z);

            #Compute the expressions for the boundary conditions of Sdot and A, to be used later.
            if fullind == 1 && EqNum == 2
                srcPrescribed = (DS0(t)/36)*(18*M*(p2-M*X^2)+18*a4-5*M^4)+(M^2/144)*(32*X*DS1(t)-61*DS2(t))+3*M^2*DS1(t)^2/(16*DS0(t));
                # srcPrescribed = ((18*a4 - 5*Power(M,4) + 18*M*p2 - 18*Power(M,2)*Power(X,2))/36.)
                continue
            end

            if fullind == 1 && EqNum == 4
                srcPrescribed = a4;
                continue
            end

            if degree == 2
                F1Mat[ptind,ptind] = co1_function(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4);
            end;

            F0Mat[ptind,ptind] = co0_function(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4);
            SRCVec[ptind] = - src_function(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, t, X, p2, a4);
        end

        if degree == 2
            mat = DiffMats2[domind,:,:] + F1Mat*DiffMats[domind,:,:] + F0Mat;
        else
            mat = DiffMats[domind,:,:] + F0Mat;
        end


        #Now impose boundary conditions on the domains - in general, we will want continuity matching on the LHS
        #The homogeneous solution is linearly independent from the particular integral
        #For 2nd order ODEs we need two independent homogeneous solutions - second one ensures we can also have derivative continuity!

        

        if domind > 1 || !regBC[EqNum]

            if degree == 2 #For 2nd order ODEs we'll generically need another independent homogeneous solution
                mat2 = copy(mat);
                SRCVec2 = zeros(T, Npts);
                SRCVec2[Npts] = 1;
                mat2[Npts,1:Npts-1] = zeros(T,Npts-1);
                mat2[Npts,Npts] = 1;
                solHomogeneous2[domind,:] = mat2 \ SRCVec2;
            end

            mat[1,2:Npts] = zeros(T, Npts-1);
            mat[1,1] = 1;
            SRCVec[1] = 0; 
            solHomogeneous[domind,:] = mat \ HomoSRCVec;
            solParticular[domind,:] = mat \ SRCVec;


        else #First domain, and a regular singular point at z = 0: we only need one solution here.
            solHomogeneous[domind,:] = zeros(T, Npts);
            solParticular[domind,:] = mat \ SRCVec;
        end
    end;

    solParticular = Unroll(solParticular);
    solHomogeneous = Unroll(solHomogeneous);
    if degree == 2
        solHomogeneous2 = Unroll(solHomogeneous2);
    end

    #Domain matching. We have 4 possibilities: degree = 1 or 2, regular singular point = true or false. Annoyingly the 4 equations in the system satiate these.

    if degree == 1 #In the nth subdomain, the solution is given as y_particular + c_n * y_homogeneous
        coeffMat = zeros(T,Ndom,Ndom);
        coeffSrc = zeros(T,Ndom);
    
        if !regBC[EqNum] #Fulfilled by the 2nd equation.
            coeffMat[1,1] = solHomogeneous[1];
            coeffSrc[1] = - solParticular[1] + srcPrescribed;
        elseif regBC[EqNum] #Fulfilled by the 3rd equation. Effectively amounts to setting c0 to zero.
            coeffMat[1,1] =  T(1.);
            coeffSrc[1] = T(0.);
        end

        for ii in 2:Ndom
            coeffMat[ii,ii - 1] = solHomogeneous[(ii - 1)*Npts];
            coeffMat[ii,ii] = -solHomogeneous[(ii-1)*Npts + 1];
            coeffSrc[ii] = solParticular[(ii-1)*Npts + 1]-solParticular[(ii-1)*Npts];
        end

        coeffVec = coeffMat \ coeffSrc;
        for ii in 1:Ndom
            locrange = (ii - 1)*Npts + 1 : ii*Npts;
            sol[locrange] = solParticular[locrange] + coeffVec[ii] .* solHomogeneous[locrange];
        end

    elseif degree == 2 #In the nth subdomain, the solution is given as y_particular + c_n * y_homogeneous + d_n * y_homogeneous2
        #Store first all c's, then all d's. Similarly first all continuity equations, then all derivative continuity equations.
        #We also need to compute the derivatives of the basis solutions.

        coeffMat = zeros(T,2*Ndom,2*Ndom);
        coeffSrc = zeros(T,2*Ndom);

        DsolParticular = ComputeSingleDerivative(solParticular);
        DsolHomogeneous = ComputeSingleDerivative(solHomogeneous);
        DsolHomogeneous2 = ComputeSingleDerivative(solHomogeneous2);

        #We may always set d0 = 0

        coeffMat[Ndom+1,Ndom+1] = T(1.);
        coeffSrc[Ndom+1] = T(0.);

        if !regBC[EqNum] #Satisfied by the 4th equation. c0 is necessary.
            coeffMat[1,1] = solHomogeneous[1];
            coeffSrc[1] = - solParticular[1] + srcPrescribed;
        elseif regBC[EqNum] #Satisfied by the 1st equation. We may also drop c0.
            coeffMat[1,1] =  T(1.);
            coeffSrc[1] = T(0.);
        end

        for ii in 2:Ndom
            #Continuity of the function:
            coeffMat[ii,ii - 1] = solHomogeneous[(ii - 1)*Npts];
            coeffMat[ii,ii] = -solHomogeneous[(ii-1)*Npts + 1];
            coeffMat[ii, Ndom + ii - 1] = solHomogeneous2[(ii - 1)*Npts];
            coeffMat[ii, Ndom + ii] = -solHomogeneous2[(ii-1)*Npts + 1];
            coeffSrc[ii] = solParticular[(ii-1)*Npts + 1]-solParticular[(ii-1)*Npts];

            #Continuity of the derivative:
            coeffMat[Ndom + ii,ii - 1] = DsolHomogeneous[(ii - 1)*Npts];
            coeffMat[Ndom + ii,ii] = -DsolHomogeneous[(ii-1)*Npts + 1];
            coeffMat[Ndom + ii, Ndom + ii - 1] = DsolHomogeneous2[(ii - 1)*Npts];
            coeffMat[Ndom + ii, Ndom + ii] = -DsolHomogeneous2[(ii-1)*Npts + 1];
            coeffSrc[Ndom + ii] = DsolParticular[(ii-1)*Npts + 1]-DsolParticular[(ii-1)*Npts];
        end

        coeffVec = coeffMat \ coeffSrc;


        for ii in 1:Ndom #Assemble the full solution, with the matching coefficients.
            locrange = (ii - 1)*Npts + 1 : ii*Npts;
            sol[locrange] = solParticular[locrange] + coeffVec[ii] .* solHomogeneous[locrange] + coeffVec[Ndom+ii] .* solHomogeneous2[locrange];
        end

    end

    return sol
end

function UnSubSdot(SdotSub, X, z, t) 
    #Computes Sdot given SdotTilde at a point
    LN = -log(z);

    sd0 = DS0(t)/2;
    sd1 = DS0(t)*X + DS1(t);
    sd2 = - (M^2 * DS0(t)^2 - 3 * DS0(t)^2 * X^2 - 6 * X * DS0(t) * DS1(t) - 3 * DS1(t)^2) / (6 * DS0(t));
    sd3 = - 2* M^2 * DS1(t) / 9;

    sigma4 = - (3*M^2*DS1(t)^2 - M^2*DS0(t)*DS2(t)) / (12 * DS0(t));
    sigma5 = M^2 * (5*DS1(t)^3 + 3*DS0(t)*DS1(t)*(5*X*DS1(t) - 3*DS2(t))+DS0(t)^2*(-5*X*DS2(t)+2*DS3(t))) / (30 * DS0(t)^2);

    sdot = sd0 / z^2 + sd1 / z + sd2 + z * sd3 + sigma4 * LN / z^2 + sigma5 * LN / z^3 + z^2 * SdotSub;
    return sdot
end

function CorrectXi(Var, X, t, margin::Int64 = 10)
    # One step in the procedure to fix ξ such that the apparent horizon is at zAH
    # Needs to be iterated over
    Xnew = copy(T(X));

    z = grid[margin];
    SdotSub = Var[3,margin];
    SdotFull = [UnSubSdot(SdotSub, Xnew, z, t)];
    subgrid = [z];



    for ii in (margin+1):N
        if ii%Npts == 1
            continue
        else
            z = grid[ii];
            SdotSub = Var[3, ii];
            push!(subgrid, z);
            push!(SdotFull,UnSubSdot(SdotSub,Xnew,z,t));
        end
    end 

    interpolant = poly.fit(subgrid, SdotFull, 20);
    foo(x,p) = interpolant(x);
    # problem = nonlin.IntervalNonlinearProblem(foo, (subgrid[1], subgrid[end]));
    problem = nonlin.NonlinearProblem(foo, zAH);
    sol = nonlin.solve(problem);
    Xnew = Xnew - 1. /(2. * zAH) + 1. /(2. * sol.u);
    
    return(Xnew);
    
end;

function PlotSdot(Var, X, t,margin::Int64 = 10)

    z = grid[margin];
    SdotSub = Var[3,margin];
    SdotFull = [UnSubSdot(SdotSub, X, z,t)];
    subgrid = [z];

    for ii in (margin+1):N
        if ii%Npts == 1
            continue
        else
            z = grid[ii];
            SdotSub = Var[3, ii];
            push!(subgrid, z);
            push!(SdotFull,UnSubSdot(SdotSub,X,z,t));
        end
    end 

    interpolant = poly.fit(subgrid, SdotFull, 20);
    foo(x,p) = interpolant(x);

    fig = plt.plot(interpolant, subgrid[1], subgrid[end],label = "Interpolant");
    plt.scatter!(fig, subgrid, SdotFull, label = "Sdot");
    plt.vline!(fig, [zAH], color = :red, linestyle = :dash, label = "zAH")
    plt.hline!(fig, [0], color = :red, linestyle = :dash, label = "y = 0")

    return(fig);
end;

function ComputeBulk(PhiArr,X,a4,t)
    # PhiTilde, ξ and a_4 together determine all data on a given time slice. This just computes it.
    localVar = zeros(T, NVar, N);
    localVar[1,1:N] = copy(PhiArr);

    for ii in 1:4
        localVar[ii+1,1:N] = LinearSolveODE(localVar, ii, a4, X,t);
    end

    return localVar;
end

function TimeDer(Var,X,t, margin)
    #Compute the time derivatives of PhiTilde, ξ and a_4. 
    #There was an issue of PhiTilde 'coming apart' at the junctions between domains.
    #Now solved by the use of polynomial interpolation. The degree of the interpolant is completely arbitrary - lots of space for experimentation
    #But setting deg = N is bad.
    p2 = Var[1,1];
    a4 = Var[5,1];
    deg = 10;
    PhiT = zeros(T,N);

    #= Begin by computing ξ'(t) at z = zAH by demanding that the horizon stay at fixed z. =#
    subgrid = deleteat!(copy(grid),Npts:Npts:(N-Npts));

    PhiArr, SArr, SdotArr, PhidotArr, AArr = eachrow(Var);
    VarZ,VarZZ = ComputeDerivatives(Var);

    subPhiArr = deleteat!(copy(PhiArr)[:],Npts:Npts:(N-Npts));
    subSArr = deleteat!(copy(SArr)[:],Npts:Npts:(N-Npts));
    subSdotArr = deleteat!(copy(SdotArr)[:],Npts:Npts:(N-Npts));
    subPhidotArr = deleteat!(copy(PhidotArr)[:],Npts:Npts:(N-Npts));
    subAArr = deleteat!(copy(AArr)[:],Npts:Npts:(N-Npts));

    Phifun = poly.fit(subgrid, subPhiArr, deg);
    Sfun = poly.fit(subgrid, subSArr, deg);
    Sdotfun = poly.fit(subgrid, subSdotArr, deg);
    Phidotfun = poly.fit(subgrid, subPhidotArr, deg);
    Afun = poly.fit(subgrid, subAArr, deg);

    PhiZfun = poly.derivative(Phifun);
    SdotZfun = poly.derivative(Sdotfun);
    AZfun = poly.derivative(Afun);

    # S = SArr[end];
    # Sdot = SdotArr[end];
    # SdotZ = VarZ[3,end];
    # Phidot = PhidotArr[end];
    # A = AArr[end];
    # AZ = VarZ[5,end];

    S = Sfun(zAH);
    Sdot = Sdotfun(zAH);
    SdotZ = SdotZfun(zAH);
    Phidot = Phidotfun(zAH);
    A = Afun(zAH);
    AZ = AZfun(zAH);

    XPrime = DtX(S, Sdot,SdotZ, Phidot, A, AZ, zAH, -log(zAH), X, p2, t);


    #= Now compute ∂_t Φ using the definition of Φdot =#


    for ii in 2:N
        z = grid[ii];
        LN = -log(z);

        # if ii > margin && ii < N - margin
        #     Phi = Phifun(z);
        #     PhiZ = PhiZfun(z);
        #     Phidot = Phidotfun(z);
        #     A = Afun(z);
        # else
        #     Phi = PhiArr[ii];
        #     PhiZ = VarZ[1,ii];
        #     Phidot = PhidotArr[ii];
        #     A = AArr[ii];
        # end

        Phi = Phifun(z);
        PhiZ = PhiZfun(z);
        Phidot = Phidotfun(z);
        A = Afun(z);

        # Phi = PhiArr[ii];
        # PhiZ = VarZ[1,ii];
        # Phidot = PhidotArr[ii];
        # A = AArr[ii];

        PhiT[ii] = DtPhi(Phi, PhiZ, Phidot, A, z, LN, X, XPrime, t);
    end

    #Use a polynomial fit to extend to the origin, otherwise get bad behaviour - place for a fix?

    tempfun = poly.fit(grid[2:15],PhiT[2:15],10);
    PhiT[1] = tempfun(grid[1]);

    # Finally compute a4'(t) using the explicit formula  

    a4Prime = Dta4(X, a4, p2, XPrime, PhiT[1], t);
    return XPrime, PhiT, a4Prime;
end;

#DIFFERENT TIME INTEGRATORS:
#'margin' is not currently used, I will get rid of it at some point.

function RK4(Var,X,a4, t,dt, margin)
    Var0 = copy(Var);
    Phi0 = Var0[1,1:N];
    X0 = copy(X);
    a40 = copy(a4);

    k1X, k1Phi, k1a4 = TimeDer(Var0, X0,t,margin);

    X1 = X0 + dt * k1X /2;
    Phi1 = Phi0 + dt * k1Phi /2;
    a41 = a40 + dt * k1a4 /2;

    Var1 = ComputeBulk(Phi1, X1, a41,t + dt/2);

    k2X, k2Phi, k2a4 = TimeDer(Var1, X1, t + dt/2,margin);

    X2 = X0 + dt * k2X /2;
    Phi2 = Phi0 + dt * k2Phi /2;
    a42 = a40 + dt * k2a4 /2;

    Var2 = ComputeBulk(Phi2,X2,a42,t + dt/2);

    k3X, k3Phi, k3a4 = TimeDer(Var2, X2,t + dt/2, margin);

    X3 = X0 + dt * k3X ;
    Phi3 = Phi0 + dt * k3Phi ;
    a43 = a40 + dt * k3a4 ;    

    Var3 = ComputeBulk(Phi3,X3,a43,t+dt);

    k4X, k4Phi, k4a4 = TimeDer(Var3, X3,t + dt, margin);

    kX = dt*(k1X + 2*k2X + 2*k3X +k4X)/6;
    kPhi = dt*(k1Phi + 2*k2Phi + 2*k3Phi +k4Phi)/6;
    ka4 = dt*(k1a4 + 2*k2a4 + 2*k3a4 +k4a4)/6;

    XNew = X0 + kX; PhiNew = Phi0 + kPhi; a4New = a40 + ka4;

    VarNew = ComputeBulk(PhiNew, XNew, a4New,t+dt);

    return VarNew, XNew, a4New

end;

function RK2(Var, X, a4, t, dt, margin)
    Var0 = copy(Var)
    Phi0 = Var0[1, 1:N]
    X0 = copy(X)
    a40 = copy(a4)

    # First derivative (k1)
    k1X, k1Phi, k1a4 = TimeDer(Var0, X0, t, margin)

    # Predictor step (Euler estimate)
    X1 = X0 + dt * k1X
    Phi1 = Phi0 + dt * k1Phi
    a41 = a40 + dt * k1a4

    Var1 = ComputeBulk(Phi1, X1, a41, t + dt)

    # Second derivative (k2)
    k2X, k2Phi, k2a4 = TimeDer(Var1, X1, t + dt, margin)

    # Heun’s RK2 update: average of k1 and k2
    kX = dt * (k1X + k2X) / 2
    kPhi = dt * (k1Phi + k2Phi) / 2
    ka4 = dt * (k1a4 + k2a4) / 2

    XNew = X0 + kX
    PhiNew = Phi0 + kPhi
    a4New = a40 + ka4

    VarNew = ComputeBulk(PhiNew, XNew, a4New, t + dt)

    return VarNew, XNew, a4New
end

function ForwardEuler(Var,X,a4, t,dt, margin)
    Var0 = copy(Var);
    Phi0 = Var0[1,1:N];
    X0 = copy(X);
    a40 = copy(a4);

    k1X, k1Phi, k1a4 = TimeDer(Var0, X0,t,margin);

    XNew = X0 + k1X * dt; PhiNew = Phi0 + k1Phi * dt; a4New = a40 + k1a4 * dt;

    VarNew = ComputeBulk(PhiNew, XNew, a4New,t+dt);

    return VarNew, XNew, a4New

end;

function BackwardEuler(Var, X, a4, t, dt, margin; tol=1e-8, maxiter=20)
    Var0 = copy(Var)
    Phi0 = Var0[1, 1:N]
    X0 = copy(X)
    a40 = copy(a4)

    # Initial guess: forward Euler step
    k1X, k1Phi, k1a4 = TimeDer(Var0, X0, t, margin)
    X_guess = X0 + dt * k1X
    Phi_guess = Phi0 + dt * k1Phi
    a4_guess = a40 + dt * k1a4

    # Fixed-point iteration for implicit update
    for iter in 1:maxiter
        Var_guess = ComputeBulk(Phi_guess, X_guess, a4_guess, t + dt)
        kX, kPhi, ka4 = TimeDer(Var_guess, X_guess, t + dt, margin)

        # New estimates
        X_new = X0 + dt * kX
        Phi_new = Phi0 + dt * kPhi
        a4_new = a40 + dt * ka4

        # Convergence check
        err = maximum(abs.([X_new - X_guess; Phi_new - Phi_guess; a4_new - a4_guess]))
        if err < tol
            X_guess, Phi_guess, a4_guess = X_new, Phi_new, a4_new
            break
        end

        # Update guess for next iteration
        X_guess, Phi_guess, a4_guess = X_new, Phi_new, a4_new
    end

    # Final update
    VarNew = ComputeBulk(Phi_guess, X_guess, a4_guess, t + dt)

    return VarNew, X_guess, a4_guess
end

function AB4(Var,X,a4, t,dt, OldTimeDer)
    Var0 = copy(Var);
    Phi0 = Var0[1,1:N];
    X0 = copy(X);
    a40 = copy(a4);

    OldF = copy(OldTimeDer);

    k0X, k0Phi, k0a4 = TimeDer(Var0, X0,t,margin);

    k1X, k1Phi, k1a4 = OldF[1];
    k2X, k2Phi, k2a4 = OldF[2];
    k3X, k3Phi, k3a4 = OldF[3];

    kX = dt * (55 * k0X - 59 * k1X + 37 * k2X - 9 * k3X) / 24;
    kPhi = dt * (55 * k0Phi - 59 * k1Phi + 37 * k2Phi - 9 * k3Phi) / 24;
    ka4 = dt * (55 * k0a4 - 59 * k1a4 + 37 * k2a4 - 9 * k3a4) / 24;


    OldTimeDer[1] = [k0X, k0Phi, k0a4];
    OldTimeDer[2] = OldF[1];
    OldTimeDer[3] = OldF[2];

    XNew = X0 + kX; PhiNew = Phi0 + kPhi; a4New = a40 + ka4;

    VarNew = ComputeBulk(PhiNew, XNew, a4New,t+dt);

    return VarNew, XNew, a4New
end

function Monitor(Var, t, X)
    # alpha is set to zero
    # Expressions copied from the paper, could have typos!
    beta = T(1/16);
    p2 = Var[1,1];
    a4 = Var[5,1];

    Eps = -3*a4/4 - M*p2 + M^2 * X^2 - M^4 * (beta - T(7/36))+3*DS1(t)^4/(16*DS0(t)^4) + M^2*(DS1(t)^2 / (8*DS0(t)^2) + 2*DS2(t)/(3*DS0(t)));
    Mom = -a4/4 + M * p2/3 - M^2 * X^2 /3 +M^4 * (beta - T(5/108)) + DS1(t)^2 * (DS1(t)^2 - 4*DS0(t)DS2(t))/(16*DS0(t)^4) - (M^2/3)*(DS1(t)^2/(8*DS0(t)^2)+13*DS2(t)/(12*DS0(t)));
    Op = -2*p2 + M * X^2 - M^3 * (4*beta - T(1/3)) + M*(-DS1(t)^2 / (4*DS0(t)^2)+5*DS2(t)/(4*DS0(t)));

    return Eps, Mom, Op
end

function AdjustGauge(Var,X,t)
    locX = copy(X);
    locVar = copy(Var);
    locp2 = locVar[1,1];
    loca4 = locVar[5,1];

    err = 1;

    while err > 1.e-15;
        locVar[2,1:N] = LinearSolveODE(locVar, 1, loca4, locX, t);
        locVar[3,1:N] = LinearSolveODE(locVar, 2, loca4, locX, t);
        locXnew = CorrectXi(locVar, locX, t, 10);
        err = abs(locXnew - locX);
        locX = locXnew;
    end
    NewVar = ComputeBulk(locVar[1,:],locX,loca4,t);
    NewX = locX;

    return NewVar, NewX
end

function Evolve(initVar, initX, inita4, inittime, maxtime, dt, write_out ,out_arr, out_monitor)
    #The full time evolution function
    #Currently set up to use AB4 but can be changed

    OldTimeDer = [];
    time = inittime;

    VarOld = copy(initVar);
    XOld = copy(initX);
    a4Old = copy(inita4);

    Eps, Mom, Op = Monitor(VarOld,time, XOld);


    push!(out_arr,[inittime, XOld, a4Old, VarOld]);
    push!(out_monitor, [Eps,Mom,Op];)

    counter = 0;

    prog = ProgressUnknown(desc = "Starting the evolution", spinner = true)

    #First couple of steps to set up later AB4
    for ii in 1:3
        push!(OldTimeDer, TimeDer(VarOld, XOld, time, 10));

        VarOld, XOld, a4Old = RK4(VarOld, XOld, a4Old, time, dt, 0);

        time = time+dt;
        counter += 1; 

        testsdot = UnSubSdot(VarOld[3,end],XOld,zAH,time);
        next!(prog, desc = string("t = ",format(time, precision=3),",  Sdot at zAH = ", format(testsdot, precision=3)))

        if counter == write_out
            Eps, Mom, Op = Monitor(VarOld, time, XOld)

            push!(out_arr,[time, XOld, a4Old, VarOld]);
            push!(out_monitor, [Eps,Mom,Op];)
            counter = 0;
        end

    end

    while time<maxtime

        VarOld, XOld, a4Old = AB4(VarOld, XOld, a4Old, time, dt, OldTimeDer);

        time = time+dt;
        counter += 1; 

        testsdot = UnSubSdot(VarOld[3,end],XOld,zAH,time);

        # if abs(testsdot)>.00004
        #     VarOld, XOld = AdjustGauge(VarOld,XOld,time);
        #     testsdot = UnSubSdot(VarOld[3,end],XOld,zAH,time);
        # end

        next!(prog, desc = string("time = ",format(time, precision=3),", Sdot at zAH = ", format(testsdot, precision=5))) 

    
        if counter == write_out
            Eps, Mom, Op = Monitor(VarOld, time, XOld)

            push!(out_arr,[time, XOld, a4Old, VarOld]);
            push!(out_monitor, [Eps,Mom,Op];)
            counter = 0;
        end
    end
end


# function ComputeEquation(Var, EqNum, a4, X)

#     VarZ, VarZZ = ComputeDerivatives(Var);
#     p2 = Var[1,1];

#     eqvec = zeros(T,N);
#     for ii in 1:N
#         Phi, S, Sdot, Phidot, A = Var[1:NVar,ii];
#         PhiZ, SZ, SdotZ, PhidotZ, AZ = VarZ[1:NVar,ii]
#         PhiZZ, SZZ, SdotZZ, PhidotZZ, AZZ = VarZZ[1:NVar,ii];
#         z = grid[ii]; LN = -T(log(z));
#         if EqNum == 1
#             eq = Equation1(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, X);
#             eqvec[ii] = eq;
#         elseif EqNum == 2
#             if ii == 1
#                 eq = Sdot - ((18*a4 - 5*Power(M,4) + 18*M*p2 - 18*Power(M,2)*Power(X,2))/36.);
#             else
#                 eq = Equation2(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, X);
#             end
#             eqvec[ii] = eq;
#         elseif EqNum == 3
#             eq = Equation3(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, X);
#             eqvec[ii] = eq;
#         elseif EqNum == 4
#             if ii == 1
#                 eq = A - a4;
#             else
#                 eq = Equation4(Phi,PhiZ,PhiZZ, S,SZ,SZZ, Sdot,SdotZ,SdotZZ, Phidot,PhidotZ,PhidotZZ, A,AZ,AZZ, z,LN, X);
#             end
#             eqvec[ii] = eq;
#         end

#     end
#     return(eqvec);
# end;

# function ComputeJacobian(Var, EqNum, a4, X)

#     eps = T(1e-8);
#     Jac = zeros(T,N,N); 

#     for ii in 1:N
#         tmpVar = copy(Var);
#         tmpVar[EqNum + 1, ii] += eps;
#         EqP = ComputeEquation(tmpVar, EqNum, a4, X);

#         tmpVar = copy(Var);
#         tmpVar[EqNum + 1, ii] -= eps;
#         EqM = ComputeEquation(tmpVar, EqNum, a4, X);

#         for kk in 1:N
#             Jac[kk,ii] = (EqP[kk] - EqM[kk])/(2*eps);
#         end
#     end
#     return(Jac);
# end;

# function NewtonSolver(Var, EqNum, a4, X)
#     LocVar = copy(Var);
#     iter = 1;
#     eqvec = ComputeEquation(LocVar, EqNum, a4, X);
#     err = linalg.norm(eqvec);
#     print("\rInitial error = ", err, " ."); 
#     flush(stdout);
#     if err < ErrorTolerance
#             flush(stdout);
#             return LocVar;
#     end

#     while iter <= IterMax
#         Jac = ComputeJacobian(LocVar, EqNum, a4, X);
#         step = Jac \ eqvec;
#         LocVar[EqNum+1, 1:N] = LocVar[EqNum+1, 1:N] - step;
#         eqvec = ComputeEquation(LocVar, EqNum, a4, X);
#         err = linalg.norm(eqvec);

#         print("\rAfter iteration n. $iter,  error = $err.");
#         flush(stdout);
#         if err < ErrorTolerance
#             flush(stdout);
#             break
#         end
#         iter += 1;

#     end
#     return LocVar;
# end;


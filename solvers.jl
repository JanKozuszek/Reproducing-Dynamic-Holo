# Functions for computing equations etc. Order of variables: Phi, S, Sdot, Phidot, A.
# X = ξ. XPrime = ξ'. In general, spatial derivatives are denoted by adding a Z at the end, i.e. AZZ = ∂_zz(A).
# The data needed at the start of a timestep: a profile for Φ (which includes ϕ_2), values for a_4 and ξ.

toFormat(x) = x             # fallback for non-collections
toFormat(x::Number) = T(x)
toFormat(x::AbstractArray) = [toFormat(y) for y in x]

# Properties of the ODEs and indices of the different variables, as needed.

global const varnamelist = ["S","Sdot","Phidot","A"]; 
global const degreelist = [2,1,1,2];
global const regBC = [true, false, true, false];
global const E = exp(one(T));

global const ind_phi    = 1;
global const ind_s      = 2;
global const ind_sdot   = 3;
global const ind_phidot = 4;
global const ind_a      = 5;

global const ind_t  = 1;
global const ind_X  = 2;
global const ind_p2 = 3;
global const ind_a4 = 4;
global const ind_S0 = 5;
global const ind_S1 = 6;

global const NVar = 5;
global const DConst = 0;

Power(x,y) = x^y;

# --------------------------------------------------------------
#  Dictionary for accessing different coefficient functions
# --------------------------------------------------------------

function_map = Dict(
    "S" => Dict(
        :Coeff0 => SCoeff0,
        :Coeff1 => SCoeff1,
        :Coeff2 => SCoeff2,
        :Src    => SSrc
    ),
    "Sdot" => Dict(
        :Coeff0 => SdotCoeff0,
        :Coeff1 => SdotCoeff1,
        :Coeff2 => SdotCoeff2,
        :Src    => SdotSrc
    ),
    "Phidot" => Dict(
        :Coeff0 => PhidotCoeff0,
        :Coeff1 => PhidotCoeff1,
        :Coeff2 => PhidotCoeff2,
        :Src    => PhidotSrc
    ),
    "A" => Dict(
        :Coeff0 => ACoeff0,
        :Coeff1 => ACoeff1,
        :Coeff2 => ACoeff2,
        :Src    => ASrc
    )
);

import LinearAlgebra as linalg
import Plots as plt
import Polynomials as poly
import NonlinearSolve as nonlin
import LsqFit as lsq

using Format
using .Threads
using ProgressMeter
using Parameters: @unpack
using StaticArrays
using FFTW
using GenericLinearAlgebra

function ComputeDerivatives(Var)
    #Takes a full state vector of dimensions NVar × N, returns all first and second z-derivatives.
    VarZ  = [similar(v) for v in Var]
    VarZZ = [similar(v) for v in Var]

    
    Threads.@threads for dom in 1:Ndom
        DM1 = diff1_mats[dom];
        DM2 = diff2_mats[dom];
        for var in 1:NVar
            VarZ[var][dom]  = DM1*Var[var][dom];
            VarZZ[var][dom] = DM2*Var[var][dom];
        end
    end

    return VarZ, VarZZ
end

function ComputeSingleDerivative(Vec; deg = 1)
    #Takes a single variable vector, returns its first, second or 'damping' derivative only.
    DVec = [similar(v) for v in Vec];

    Threads.@threads for dom in 1:Ndom
        DM = (deg <=2 ? (deg ==1 ? diff1_mats[dom] : diff2_mats[dom]) : damp_mats[dom]);
        DVec[dom] = DM*Vec[dom];
    end
    return DVec
end

function BoundaryInterpolate(VarVec)
    # Take the near - boundary points in the given variable and fit a model, including some logarithmic terms, to extract the FG expansion coefficients
    model(z, p) = p[1] .+ p[2] .*z + p[3] .*z.^2 .+ p[4]  .* z.^3 .+ p[5] .* log.(z .+eps(T)) .* z .+ p[6] .* log.(z .+eps(T)) .* z .^2;

    zdata = grids[1][1:7];
    ydata = VarVec[1:7];

    p0 = [VarVec[1], 0., 0., 0. , 0., 0.];
    fit = lsq.curve_fit(model, zdata, ydata, p0);
    params = lsq.coef(fit);
    return params
end

#---------------------------------------------
#ODE SOLVER FUNCTIONS
#---------------------------------------------

function ComputeODEMatrix(EqNum, domind, params, Var, VarZ, VarZZ, p3,p4)
    #Compute the matrix for one of the ODEs, so that it can be written as mat * vec = src.
    #Working in a given domain "domind"
    #Not dealing with boundary conditions here.

    degree = degreelist[EqNum];
    var = varnamelist[EqNum];

    t, X, p2, a4, DS0, DS1 = params;
    Npts = domain_sizes[domind];

    # srcPrescribed = zero(T);

    if degree == 2
        co2_function = function_map[var][:Coeff2]
    end

    co1_function = function_map[var][:Coeff1]
    co0_function = function_map[var][:Coeff0]
    src_function = function_map[var][:Src]


    if degree == 2
        F2Vec = Vector{T}(undef, Npts);
    end
    F1Vec = Vector{T}(undef, Npts);
    F0Vec = Vector{T}(undef, Npts);
    SRCVec = Vector{T}(undef, Npts);

    mat = Matrix{T}(undef,Npts,Npts)

    DS2 = DS2f(params, p3, p4);
    DS3 = DS3f(params, p3, p4);
    DS4 = DS4f(params, p3, p4);

    for ptind in 1:Npts
        vals   =  [x[domind][ptind] for x in Var];
        valsZ  =  [x[domind][ptind] for x in VarZ];
        valsZZ =  [x[domind][ptind] for x in VarZZ];

        z = grids[domind][ptind]
        LN = ifelse(z == 0, -log(eps(typeof(z))), -log(z)) #Prevents evaluating log(0);

        pdevars = PDEVars(vals..., valsZ..., valsZZ..., z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

        if degree == 2
            F2Vec[ptind] = co2_function(pdevars);
        end

        F1Vec[ptind] = co1_function(pdevars);
        F0Vec[ptind] = co0_function(pdevars);
        SRCVec[ptind] = - src_function(pdevars);

    end

    DM1 = diff1_mats[domind];   
    DM2 = (degree==2 ?  diff2_mats[domind] : nothing);


    if degree == 2
        mat = F2Vec .* DM2 + F1Vec .* DM1 + linalg.Diagonal(F0Vec);
    else
        mat = F1Vec .* DM1  + linalg.Diagonal(F0Vec);
    end

    # if isnan(SRCVec[1])||isinf(SRCVec[1])
    #     SRCVec[1] = BoundaryInterpolate(SRCVec)[1];
    # end
    return mat, SRCVec

end

function DomainMatching(EqNum, degree, solParticular, solHomogeneous, solHomogeneous2)
    #Domain matching. We have 4 possibilities: degree = 1 or 2, regular singular point = true or false. Annoyingly the 4 equations in the system satiate these.
    #Now moved dealing with the BCs to the matrix expressions themselves.

    sol = deepcopy(zero_var);

    if degree == 1 #In the nth subdomain, the solution is given as y_particular + c_n * y_homogeneous
        coeffMat = zeros(T,Ndom,Ndom);
        coeffSrc = zeros(T,Ndom);

        coeffMat[1,1] =  T(1.);
        coeffSrc[1] = T(0.);

        for dom in 2:Ndom
            coeffMat[dom,dom - 1] = solHomogeneous[dom-1][end];
            coeffMat[dom,dom] = -solHomogeneous[dom][1];
            coeffSrc[dom] = solParticular[dom][1]-solParticular[dom-1][end];
        end

        coeffVec = coeffMat \ coeffSrc;

        for dom in 1:Ndom
            sol[dom] = solParticular[dom] + coeffVec[dom] .* solHomogeneous[dom];
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

        coeffMat[1,1] =  T(1.);
        coeffSrc[1] = T(0.);

        for dom in 2:Ndom
            coeffMat[dom,dom - 1] = solHomogeneous[dom-1][end];
            coeffMat[dom,dom] = -solHomogeneous[dom][1];
            coeffSrc[dom] = solParticular[dom][1]-solParticular[dom-1][end];


            #Continuity of the function:
            coeffMat[dom,dom - 1] = solHomogeneous[dom - 1][end];
            coeffMat[dom,dom] = -solHomogeneous[dom][1];
            coeffMat[dom, Ndom + dom - 1] = solHomogeneous2[dom-1][end];
            coeffMat[dom, Ndom + dom] = -solHomogeneous2[dom][1];
            coeffSrc[dom] = solParticular[dom][1]-solParticular[dom-1][end];

            #Continuity of the derivative:
            coeffMat[Ndom + dom,dom - 1] = DsolHomogeneous[dom-1][end];
            coeffMat[Ndom + dom,dom] = -DsolHomogeneous[dom][1];
            coeffMat[Ndom + dom, Ndom + dom - 1] = DsolHomogeneous2[dom-1][end];
            coeffMat[Ndom + dom, Ndom + dom] = -DsolHomogeneous2[dom][1];
            coeffSrc[Ndom + dom] = DsolParticular[dom][1]-DsolParticular[dom-1][end];
        end

        coeffVec = coeffMat \ coeffSrc;

        for dom in 1:Ndom #Assemble the full solution, with the matching coefficients.
            sol[dom] = solParticular[dom] + coeffVec[dom] .* solHomogeneous[dom] + coeffVec[Ndom+dom] .* solHomogeneous2[dom];
        end

    end

    return sol
end

function LinearSolveODE(EqNum, params, Var, VarZ, VarZZ)
    # Solve an ODE by direct linear inversion. Grid is subdivided into Ndom domains of Npts each.
    # This routine can now handle the case where zmin = 0. But that might still be bad for stability.
    degree = degreelist[EqNum];

    sol = deepcopy(zero_var);
    solParticular = deepcopy(zero_var);
    solHomogeneous = deepcopy(zero_var);
    if degree == 2 
        solHomogeneous2 = deepcopy(zero_var);
    end

    p3, p4 = BoundaryInterpolate(Var[1][1])[2:3];

    # srcPrescribed = zero(T);

    #Basic parallelization here:
    Threads.@threads for domind in 1:Ndom
        Npts = domain_sizes[domind];

        
        mat, SRCVec = ComputeODEMatrix(EqNum, domind, params, Var, VarZ, VarZZ, p3, p4);

        HomoSRCVec = zeros(T, Npts);
        HomoSRCVec[1] = 1;
    
        
        #Now impose boundary conditions on the domains - in general, we will want continuity matching on the LHS
        #The homogeneous solution is linearly independent from the particular integral
        #For 2nd order ODEs we need two independent homogeneous solutions - second one ensures we can also have derivative continuity!        

        if domind > 1

            if degree == 2 #For 2nd order ODEs we'll generically need another independent homogeneous solution
                mat2 = copy(mat);
                SRCVec2 = zeros(T, Npts);
                SRCVec2[Npts] = 1;
                mat2[Npts,1:Npts-1] = zeros(T,Npts-1);
                mat2[Npts,Npts] = 1;
                solHomogeneous2[domind] = mat2 \ SRCVec2;
            end

            SRCVec[1] = 0; 
            mat[1,2:Npts] = zeros(T, Npts-1);
            mat[1,1] = 1;

            lumat = linalg.lu(mat);
            solHomogeneous[domind] = lumat \ HomoSRCVec;
            solParticular[domind] = lumat \ SRCVec;


        else #First domain: we only need one solution here.
            # SRCVec[1] = BoundaryInterpolate(SRCVec)[1];
            solHomogeneous[domind] = zeros(T, Npts);
            solParticular[domind] = mat \ SRCVec;
        end
    end;

    sol = (degree == 1 ? DomainMatching(EqNum,degree,solParticular,solHomogeneous,nothing) : DomainMatching(EqNum,degree,solParticular,solHomogeneous,solHomogeneous2))

    return sol
end

function ComputeBulk(params, Vararg)
    # PhiTilde, ξ and a_4 together determine all data on a given time slice. This just computes it.
    Var = deepcopy(Vararg);
    VarZ, VarZZ = ComputeDerivatives(Var);

    for eq in 1:4
        Var[eq+1]= LinearSolveODE(eq, params, Var, VarZ, VarZZ);
        VarZ[eq+1] = ComputeSingleDerivative(Var[eq+1]);
        VarZZ[eq+1] = ComputeSingleDerivative(Var[eq+1], deg=2);
    end

    return Var;
end

function ComputeBulkFromVec(params, phi_vec)
    # PhiTilde, ξ and a_4 together determine all data on a given time slice. This just computes it.
    # Takes only a vector for phi as an argument!
    Var = [copy(zero_var) for var in 1:NVar];
    Var[1] = copy(phi_vec);
    VarZ, VarZZ = ComputeDerivatives(Var);

    for eq in 1:4
        Var[eq+1]= LinearSolveODE(eq, params, Var, VarZ, VarZZ);
        VarZ[eq+1] = ComputeSingleDerivative(Var[eq+1]);
        VarZZ[eq+1] = ComputeSingleDerivative(Var[eq+1], deg=2);
    end

    return Var;
end

function CorrectXi(params, Var)
    # One step in the procedure to fix ξ such that the apparent horizon is at zAH
    # Needs to be iterated over
    Xnew = copy(params[ind_X]);
    DS0 = copy(params[ind_S0]);
    DS1 = copy(params[ind_S1])

    p3, p4 = BoundaryInterpolate(Var[1][1])[2:3];
    DSVals = [DS0, DS1, DS2f(params, p3, p4), DS3f(params, p3, p4), DS4f(params, p3, p4)];

    flatgrid = reduce(vcat,[subgrid[1:end-1] for subgrid in grids[2:end]]);
    flatVar = reduce(vcat,[subvar[1:end-1] for subvar in Var[ind_sdot][2:end]]);
    Sdot_unsub = [Sdot_to_unsub(params, flatgrid[ii], -log(flatgrid[ii]), DSVals, flatVar[ii]) for ii in eachindex(flatgrid)];


    # This interpolation degree is arbitrary!!
    interpolant = poly.fit(flatgrid, Sdot_unsub, 20);
    foo(x,p) = interpolant(x);
    problem = nonlin.NonlinearProblem(foo, zAH);
    sol = nonlin.solve(problem);
    Xnew = Xnew - 1. /(2. * zAH) + 1. /(2. * sol.u);
    
    return(Xnew);
    
end;

function PlotSdot(params, Var)

    p3, p4 = BoundaryInterpolate(Var[1][1])[2:3];
    DS0 = copy(params[ind_S0]);
    DS1 = copy(params[ind_S1]);
    DSVals = [DS0, DS1, DS2f(params, p3, p4), DS3f(params, p3, p4), DS4f(params, p3, p4)];
    
    flatgrid = reduce(vcat,[subgrid[1:end-1] for subgrid in grids[2:end]]);
    flatVar = reduce(vcat,[subvar[1:end-1] for subvar in Var[ind_sdot][2:end]]);
    Sdot_unsub = [Sdot_to_unsub(params, flatgrid[ii], -log(flatgrid[ii]), DSVals, flatVar[ii]) for ii in eachindex(flatgrid)];

    interpolant = poly.fit(flatgrid, Sdot_unsub, 20);
    foo(x,p) = interpolant(x);

    fig = plt.plot(interpolant, flatgrid[1], flatgrid[end],label = "Interpolant");
    plt.scatter!(fig, flatgrid, Sdot_unsub, label = "Sdot");
    plt.vline!(fig, [zAH], color = :red, linestyle = :dash, label = "zAH")
    plt.hline!(fig, [0], color = :red, linestyle = :dot, label = "y = 0")

    return(fig);
end;


#---------------------------------------------
# TIME EVOLUTION FUNCTIONS
#---------------------------------------------

function TimeDer(params, Var)
    #Compute the time derivatives of PhiTilde, ξ and a_4. 
    #There was an issue of PhiTilde 'coming apart' at the junctions between domains.
    #Now solved by the correct treatment of junction points, as detailed in the Jecco paper. 
    t, X, p2, a4, DS0, DS1= params;
    p3, p4 = BoundaryInterpolate(Var[1][1])[2:3];

    PhiT = deepcopy(zero_var);

    DS2 = DS2f(params, p3, p4);
    DS3 = DS3f(params, p3, p4);
    DS4 = DS4f(params, p3, p4);

    #= Begin by computing ξ'(t) at z = zAH by demanding that the horizon stay at fixed z. =#
    AArr = Var[ind_a];

    VarZ,VarZZ = ComputeDerivatives(Var);

    valsAH   =  [x[domAH][indAH] for x in Var];
    valsZAH  =  [x[domAH][indAH] for x in VarZ];
    valsZZAH =  [x[domAH][indAH] for x in VarZZ];
    pdevars = PDEVars(valsAH..., valsZAH..., valsZZAH..., zAH, -log(zAH), t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

    XPrime = DtX(pdevars);


    #= Now compute ∂_t Φ using the definition of Φdot =#


    Threads.@threads for domind in 1:Ndom
        for ptind in 1:domain_sizes[domind]

            vals   =  [x[domind][ptind] for x in Var];
            valsZ  =  [x[domind][ptind] for x in VarZ];
            valsZZ =  [x[domind][ptind] for x in VarZZ];

            z = grids[domind][ptind]
            LN = ifelse(z == 0, -log(eps(typeof(z))), -log(z)) #Prevents evaluating log(0);

            pdevars = PDEVars(vals..., valsZ..., valsZZ..., z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

            PhiT[domind][ptind] = DtPhi(XPrime, pdevars);
        end

        PhiT[domind] = cheb_filter(PhiT[domind]);
    end


    # Including log terms in near boundary interpolation. 

    # PhiT[1][1] = BoundaryInterpolate(PhiT[1])[1];

    # Propagation matching
    for dom in 1:Ndom-1
        z = grids[dom][end];
        LN = -log(z);
        csign = A_to_unsub(params, z, LN, [DS0, DS1, DS2, DS3, DS4], AArr[dom][end], XPrime);
        if csign>0
            PhiT[dom][end] = PhiT[dom+1][1];
        else
            PhiT[dom+1][1] = PhiT[dom][end];
        end
    end
    # Finally compute a4'(t) using the explicit formula  

    a4Prime = Dta4(params, [DS0, DS1, DS2, DS3, DS4], XPrime, PhiT[1][1]);

    DtDS1 = DtDS1f(params, XPrime, PhiT[1][1], a4Prime);

    return XPrime, PhiT, a4Prime, DS1, DtDS1;
end;

function RK4(params, Var, dt)
    t, X, p2, a4, DS0, DS1 = params;
    
    Var0 = copy(Var);
    Phi0 = Var0[1];
    X0 = copy(X);
    a40 = copy(a4);

    k1X, k1Phi, k1a4, k1S0, k1S1 = TimeDer([t, X0, Phi0[1][1], a40, DS0, DS1],Var0);

    X1 = X0 + dt * k1X /2;
    Phi1 = Phi0 + dt * k1Phi /2;
    a41 = a40 + dt * k1a4 /2;
    DS01 = DS0 + dt * k1S0/2;
    DS11 = DS1 + dt * k1S1/2;

    Var1 = ComputeBulkFromVec([t+dt/2, X1, Phi1[1][1], a41, DS01, DS11],Phi1);

    k2X, k2Phi, k2a4, k2S0, k2S1 = TimeDer([t+dt/2, X1, Phi1[1][1], a41, DS01, DS11],Var1);

    X2 = X0 + dt * k2X /2;
    Phi2 = Phi0 + dt * k2Phi /2;
    a42 = a40 + dt * k2a4 /2;
    DS02 = DS0 + dt * k2S0 /2;
    DS12 = DS1 + dt * k2S1 / 2;

    Var2 = ComputeBulkFromVec([t+dt/2, X2, Phi2[1][1], a42, DS02, DS12],Phi2);

    k3X, k3Phi, k3a4, k3S0, k3S1 = TimeDer([t+dt/2, X2, Phi2[1][1], a42, DS02, DS12],Var2);

    X3 = X0 + dt * k3X ;
    Phi3 = Phi0 + dt * k3Phi ;
    a43 = a40 + dt * k3a4 ;    
    DS03 = DS0 + dt*k3S0;
    DS13 = DS1 + dt*k3S1;

    Var3 = ComputeBulkFromVec([t+dt, X3, Phi3[1][1], a43, DS03, DS13],Phi3);

    k4X, k4Phi, k4a4, k4S0, k4S1 = TimeDer([t+dt, X3, Phi3[1][1], a43, DS0, DS13],Var3);

    kX = dt*(k1X + 2*k2X + 2*k3X + k4X)/6;
    kPhi = dt*(k1Phi + 2*k2Phi + 2*k3Phi + k4Phi)/6;
    ka4 = dt*(k1a4 + 2*k2a4 + 2*k3a4 + k4a4)/6;
    kS04 = dt*(k1S0 + 2*k2S0 + 2*k3S0 + k4S0)/6;
    kS14 = dt*(k1S1 + 2*k2S1 + 2*k3S1 + k4S1)/6;

    XNew = X0 + kX; PhiNew = Phi0 + kPhi; a4New = a40 + ka4; DS0New = DS0 + kS04; DS1New = DS1 + kS14;

    p2New = PhiNew[1][1];
    param_new = [t+dt, XNew, p2New, a4New, DS0New, DS1New];

    VarNew = ComputeBulkFromVec(param_new, PhiNew);

    return param_new, VarNew

end

function AB4(params, Var ,dt, OldTimeDer)
    
    t, X0, p2, a40, DS0, DS1 = params;
    # DS1 = DS1f(params); # <- this is more stable!


    OldF = deepcopy(OldTimeDer);
    Phi0 = deepcopy(Var[1]);

    k0X, k0Phi, k0a4, k0S0, k0S1 = TimeDer(params,Var);

    k1X, k1Phi, k1a4, k1S0, k1S1 = OldF[3];
    k2X, k2Phi, k2a4, k2S0, k2S1 = OldF[2];
    k3X, k3Phi, k3a4, k3S0, k3S1 = OldF[1];

    kX = dt * (55 * k0X - 59 * k1X + 37 * k2X - 9 * k3X) / 24;
    kPhi = dt * (55 * k0Phi - 59 * k1Phi + 37 * k2Phi - 9 * k3Phi) / 24;
    ka4 = dt * (55 * k0a4 - 59 * k1a4 + 37 * k2a4 - 9 * k3a4) / 24;
    kS0 = dt * (55 * k0S0 - 59 * k1S0 + 37 * k2S0 - 9 * k3S0) / 24;
    kS1 = dt * (55 * k0S1 - 59 * k1S1 + 37 * k2S1 - 9 * k3S1) / 24;


    OldTimeDer[3] = [k0X, k0Phi, k0a4, k0S0, k0S1];
    OldTimeDer[2] = OldF[3];
    OldTimeDer[1] = OldF[2];

    XNew = X0 + kX; PhiNew = Phi0 + kPhi; a4New = a40 + ka4; DS0New = DS0 + kS0; DS1New = DS1 + kS1;
    p2New = PhiNew[1][1];


    param_new = [t+dt,XNew, p2New, a4New, DS0New, DS1New];
    VarNew = ComputeBulkFromVec(param_new,PhiNew);

    return param_new, VarNew
end

function Evolve(initparams, initVar, maxtime, dt, write_out ,out_io, monitor_io)
    #The full time evolution function
    #Currently set up to use AB4 but can be changed

    inittime, initX, initp2, inita4, initDS0, initDS1 = initparams;

    VarCurrent = deepcopy(initVar);
    XCurrent = copy(initX);
    a4Current = copy(inita4);


    OldTimeDer = []; #For the AB4 integrator, need 3 past values
    OldSdot = []; OldXarr = [];#For monitoring the constraint, need 4 past values
    time = inittime;

    CurrentParams = copy(initparams);

    Eps, Mom, Op = Monitor(CurrentParams);

    push!(out_io, VarCurrent);
    push!(monitor_io,[time,CurrentParams,Eps,Mom,Op]);


    # out_data = tofloat64(VarCurrent);
    # out_monit = vcat(Float64(inittime), Float64(XCurrent), Float64(a4Current),Float64(Eps),Float64(Mom),Float64(Op),Float64(0.));
    # write(out_io,out_data);
    # write(monitor_io, out_monit);

    push!(OldSdot,VarCurrent[ind_sdot])
    push!(OldXarr,XCurrent);

    counter = 0;

    prog = ProgressUnknown(desc = "Starting the evolution", spinner = true)

    # First couple of steps to set up later AB4
    for ii in 1:3

        CurrentParams, VarCurrent = RK4(CurrentParams, VarCurrent, dt);

        time = time+dt;
        counter += 1; 

        push!(OldTimeDer, TimeDer(CurrentParams,VarCurrent));
        push!(OldSdot,VarCurrent[ind_sdot])
        push!(OldXarr,XCurrent);
        if counter == write_out
            push!(out_io, VarCurrent);
            push!(monitor_io,[time,CurrentParams,Eps,Mom,Op]);
            counter = 0;
        end


        # testsdot = UnSubSdot(VarOld[3,end],XOld,zAH,time);
        next!(prog, desc = string("t = ",format(time, precision=3)#=,",  Sdot at zAH = ", format(testsdot, precision=3)=#))

        # if counter == write_out
        #     Eps, Mom, Op = Monitor(CurrentParams)

        #     out_data = tofloat64(VarCurrent);
        #     out_monit = vcat(Float64(time), Float64(XCurrent), Float64(a4Current),Float64(Eps),Float64(Mom),Float64(Op),Float64(0.));
        #     write(out_io,out_data);
        #     write(monitor_io, out_monit);

        #     counter = 0;
        # end

    end

    while time<maxtime

        CurrentParams, VarCurrent = AB4(CurrentParams, VarCurrent, dt, OldTimeDer);
        # CurrentParams, VarCurrent = RK4(CurrentParams,VarCurrent, dt);

        time = time+dt;
        counter += 1; 

        constr_norm= EvaluateConstraint(CurrentParams, VarCurrent, OldSdot, OldXarr, dt);

        p3, p4 = BoundaryInterpolate(VarCurrent[1][1])[2:3];

        DS0 = CurrentParams[ind_S0]
        DS1 = CurrentParams[ind_S1]
        DS2 = DS2f(CurrentParams, p3, p4);
        DS3 = DS3f(CurrentParams, p3, p4);
        DS4 = DS4f(CurrentParams, p3, p4);

        testsdot = Sdot_to_unsub(CurrentParams, zAH, -log(zAH), [DS0, DS1, DS2, DS3, DS4], VarCurrent[ind_sdot][domAH][indAH]);

        next!(prog, desc = string("time = ",format(time, precision=3),", constraint violation = ", format(constr_norm, precision=3),",  Sdot at zAH = ", format(testsdot, precision=3)));

    
        if counter == write_out
            Eps, Mom, Op = Monitor(CurrentParams)
            push!(out_io, VarCurrent);
            push!(monitor_io,[time,CurrentParams,Eps,Mom,Op]);


            # out_data = tofloat64(VarCurrent);
            # out_monit = vcat(Float64(time), Float64(XCurrent), Float64(a4Current),Float64(Eps),Float64(Mom),Float64(Op),Float64(constr_norm));
            # write(out_io,out_data);
            # write(monitor_io, out_monit);

            # flush(out_io);
            # flush(monitor_io);

            counter = 0;
        end
    end

    return CurrentParams, VarCurrent;
end

function BackwardsTimeDerivative(a4,a3,a2,a1,a0,dt)
    res = (3*a4 - 16*a3 + 36*a2 - 48*a1 + 25*a0) / (12 * dt);
    return res
end

"""
cheb_filter(f::Vector{T})
Attempt at implementing a high-frequency filter.
"""
function cheb_filter(f::Vector{T})
    N = length(f) - 1
    alpha = 36.0437;
    # alpha = 100;
    # # Transform to Chebyshev coefficients using DCT-I
    a = dct(f, 1)

    # Construct modal filter
    n = (0:N)./N;
    sigma = exp.(-alpha * (n.^ (8*N)));

    # Apply damping
    a_filtered = a .* sigma

    # Back to physical space
    f_filtered = idct(a_filtered, 1);  # normalization for DCT-I

    return f_filtered
end


#---------------------------------------------
# CHANGING TO FULL VARIABLES AND MONITORING
#---------------------------------------------

"""
EvaluateConstraint(params, Var, previous_sdot_arr_arg, previous_x_arr_arg, dt)
"""
function EvaluateConstraint(params, Var, previous_sdot_arr_arg, previous_x_arr_arg, dt)
    previous_x_arr = copy(previous_x_arr_arg);
    previous_sdot_arr = copy(previous_sdot_arr_arg);

    t, X, p2, a4, DS0, DS1 = params;
    p3,p4 = BoundaryInterpolate(Var[1][1])[2:3];

    DS2 = DS2f(params, p3, p4);
    DS3 = DS3f(params, p3, p4);
    DS4 = DS4f(params, p3, p4);

    VarZ, VarZZ = ComputeDerivatives(Var);
    res = copy(zero_var);

    XPrime = BackwardsTimeDerivative(previous_x_arr[1],previous_x_arr[2],previous_x_arr[3],previous_x_arr[4],X,dt)

    

    for domind in 1:Ndom
        for ptind in 1:domain_sizes[domind]
            vals   =  [x[domind][ptind] for x in Var];
            valsZ  =  [x[domind][ptind] for x in VarZ];
            valsZZ =  [x[domind][ptind] for x in VarZZ];

            z = grids[domind][ptind]
            LN = ifelse(z == 0, -log(eps(typeof(z))), -log(z)) #Prevents evaluating log(0);

            pdevars = PDEVars(vals..., valsZ..., valsZZ..., z, LN, t, X, p2, a4, DS0, DS1, DS2, DS3, DS4)

            SdotT = BackwardsTimeDerivative(previous_sdot_arr[1][domind][ptind],previous_sdot_arr[2][domind][ptind],previous_sdot_arr[3][domind][ptind],previous_sdot_arr[4][domind][ptind],Var[ind_sdot][domind][ptind],dt);

            res[domind][ptind] = constr(SdotT, XPrime, pdevars);
        end
    end
    res[1][1] = 0;


    previous_sdot_arr_arg[1] = previous_sdot_arr[2];
    previous_sdot_arr_arg[2] = previous_sdot_arr[3];
    previous_sdot_arr_arg[3] = previous_sdot_arr[4];
    previous_sdot_arr_arg[4] = Var[ind_sdot];

    previous_x_arr_arg[1] = previous_x_arr[2];
    previous_x_arr_arg[2] = previous_x_arr[3];
    previous_x_arr_arg[3] = previous_x_arr[4];
    previous_x_arr_arg[4] = X;

    flatres = reduce(vcat,[subres[2:end] for subres in res]);

    return linalg.norm(flatres);
end

function Sdot_to_unsub(params, z, LN, DSVals, sub_value)

    t, X, p2, a4 = params;
    Sdot = sub_value;

    DS0, DS1, DS2, DS3, DS4 = DSVals;

    unsub_value = (30*DS1^3*LN*M^2*z^5 - 30*DS0^3*(-3 - 6*X*z + M^2*z^2 - 3*X^2*z^2) + 9*DS0*DS1*z^2*(-6*DS2*LN*M^2*z^3 + 5*DS1*(2 - LN*M^2*z^2 + 2*LN*M^2*X*z^3)) + DS0^2*z*(-20*DS1*(-9 - 9*X*z + 2*M^2*z^2) + 3*LN*M^2*z^3*(4*DS3*z + DS2*(5 - 10*X*z)) + 180*z^3*Sdot))/(180*DS0^2*z^2);
    return unsub_value

end

function A_to_unsub(params, z, LN, DSVals, sub_value, XPrime)

    t, X, p2, a4 = params;
    A = sub_value;
    DS0, DS1, DS2, DS3, DS4 = DSVals;

    unsub_value = (z^3*(DS0*(DS3*LN*M^2 + DS2*(-4*LN*M^2*X - 6/z^3 + (2*LN*M^2)/z)) + DS1*(-(DS2*LN*M^2) + (3*DS1)/z^3) + (DS0^2*(3 + 6*X*z - 2*M^2*z^2 + 3*X^2*z^2 - 6*XPrime*z^2 + 3*A*z^4))/z^5))/(3*DS0^2);

    return unsub_value

end

function Monitor(params)
    # alpha is set to zero
    # Expressions copied from the paper, could have typos!
    beta = T(1/16);
    
    t, X, p2, a4, DS0, DS1 = params;

    DS2 = DS2f(params, 0,0);

    Eps = -3*a4/4 - M*p2 + M^2 * X^2 - M^4 * (beta - T(7/36))+3*DS1^4/(16*DS0^4) + M^2*(DS1^2 / (8*DS0^2) + 2*DS2/(3*DS0));
    Mom = -a4/4 + M * p2/3 - M^2 * X^2 /3 +M^4 * (beta - T(5/108)) + DS1^2 * (DS1^2 - 4*DS0*DS2)/(16*DS0^4) - (M^2/3)*(DS1^2/(8*DS0^2)+13*DS2/(12*DS0));
    Op = -2*p2 + M * X^2 - M^3 * (4*beta - T(1/3)) + M*(-DS1^2 / (4*DS0^2)+5*DS2/(4*DS0));

    return Eps, Mom, Op
end

"""
Temperature(params, subA_arr)
Computing the temperature of the *apparent* horizon
"""
function Temperature(params, subA_arr, near_boundary_phi, XPrime)
# Only need the values for A in the subdomain in which the apparent horizon is located
    t, X, p2, a4, DS0, DS1 = params;
    p3, p4 = BoundaryInterpolate(near_boundary_phi)[2:3];

    DSVals = [DS0, DS1, DS2f(params,p3,p4), DS3f(params, p3, p4), DS4f(params, p3, p4)];

    A_arr = [A_to_unsub(params, grids[domAH][ind], -log(grids[domAH][ind]), DSVals, subA_arr[ind], XPrime) for ind in eachindex(subA_arr)];
    AZ_arr = diff1_mats[domAH]*A_arr;

    return -zAH^2 * AZ_arr[indAH]
end
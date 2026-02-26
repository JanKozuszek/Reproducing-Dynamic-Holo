using LinearAlgebra

function chebgrid(start, stop, length)
    # Compute a Chebyshev-spaced grid for arbitrary endpoints
    # x - grid between -1 and 1. y - desired grid. y = alpha*x + beta.

    start = T(start);
    stop  = T(stop);


    # Chebyshev grid on [-1, 1]
    xarr = [cos(T(pi) * i / (T(length) - T(1.))) for i in 0:Int64(length)-1]

    alpha = (stop - start) / T(2.);
    beta  = (stop + start) / T(2.);

    yarr = [alpha * x + beta for x in xarr];
    reverse!(yarr);

    return yarr

end;

function chebmat(grid::AbstractVector{T})
    # Compute a differentiation matrix given a (Chebyshev) grid
    # Formula copied from the Trefethen book

    N = length(grid);
    Dmat = zeros(T, N, N);
    for ii in 1:N
        for jj in 1:N
            ci = (ii == 1 || ii == N ? T(2.) : T(1.));
            cj = (jj == 1 || jj == N ? T(2.) : T(1.));
            if ii != jj
                Dmat[ii,jj] = T((ci/cj) * (-1.)^(ii+jj) / (grid[ii] - grid[jj]));
            end
        end
    end

    for ii in 1:N
        for jj in 1:N
            if ii != jj
                Dmat[ii,ii] -= Dmat[ii,jj];
            end
        end
    end

    return Dmat;
end;

function cheb(start, stop, length)
    grid = chebgrid(start, stop, length);
    mat = chebmat(grid);

    return([mat, grid]);
end;

function finitegrid(start, stop, length)
    return range(T(start), T(stop), length)
end

function fd_weights_first(offsets)
    n = length(offsets)
    A = zeros(T, n, n)
    b = zeros(T, n)

    # Vandermonde system enforcing polynomial exactness
    for k in 0:n-1
        A[k+1, :] .= offsets .^ k
    end

    b[2] = 1.0   # first derivative of x at 0 is 1

    return A \ b
end


function finitemat(grid::AbstractVector{T}; order = 4)
    @assert iseven(order) "order must be even"

    r = order ÷ 2; 
    N = length(grid);
    Dmat = zeros(T, N,N);
    h = grid[2] - grid[1]; #Assuming uniform grid.


    for i in 1:N
        if i <= r
            #Near left boundary
            offsets = -(i-1):(order - (i - 1))
            cols = 1:(order+1);
        elseif i > N - r
            offsets = -(order - (N - i)):(N - i)
            cols = (N - order):N

        else
            offsets = -r:r
            cols = (i-r):(i+r);
        end

        w = fd_weights_first(collect(offsets)) ./ h
        Dmat[i, cols] .= w

    end

    return Dmat
end

function finitediff(start, stop, length, OrderDiff)

    grid =finitegrid(start, stop, length);
    mat = finitemat(grid; order = OrderDiff)

    return([mat, grid])

end


function MultiGridChebyshev(start, stop, Ndom, Npts)
    #Set up a multidomain grid, all Ndom domains have the same number of points Npts and the same length - could generalize that.
    
    allmats = [];
    allmats2 = [];
    grids = [];
    dampmats = [];
    L = (stop - start)/Ndom;
    loc_start = copy(start);
    loc_stop = loc_start + L;

    for ii in 1:Ndom
        D1, grid = cheb(loc_start, loc_stop, Npts);
        push!(grids, grid);
        push!(allmats, D1);
        push!(allmats2, D1*D1);
        push!(dampmats,D1*D1*D1*D1);
        loc_start = loc_stop;
        loc_stop = loc_start+L;
    end;

    return allmats2, allmats, grids, dampmats
end

function MultiGridChebyshev(endpoints::Vector{T}, lengths::Vector{Int64})
    Ndom = length(endpoints)-1;
    if !(Ndom == length(lengths))
        println("Mismatched arrays!!");
        return 1
    end

    allmats = [];
    allmats2 = [];
    dampmats = [];
    grids = [];

    for ii in 1:Ndom
        loc_start = endpoints[ii];
        loc_stop = endpoints[ii+1];
        Npts = lengths[ii];
        DM1, grid1 = cheb(loc_start, loc_stop, Npts)
        DM2 = DM1*DM1;

        push!(allmats, DM1);
        push!(allmats2, DM2);
        push!(dampmats, DM2*DM2);
        push!(grids, grid1);

    end;

    return allmats2, allmats, grids, dampmats
end;

function MultiGridMultiMethod(endpoints::Vector{T}, lengths::Vector{Int64}, orders::Vector{Int64})
    Ndom = length(endpoints) - 1;
    @assert (Ndom == length(lengths)) "Mismatched arrays!"

    allmats = [];
    allmats2 = [];
    dampmats = [];
    grids = [];

    for ii in 1:Ndom
        loc_start = endpoints[ii];
        loc_stop = endpoints[ii+1];
        Npts = lengths[ii];

        DM1 = zeros(T, Npts, Npts);
        grid1 = zeros(T, Npts);

        if orders[ii] == 0
            DM1, grid1 = cheb(loc_start, loc_stop, Npts)
        else
            DM1, grid1 = finitediff(loc_start, loc_stop, Npts, orders[ii])
        end

        DM2 = DM1*DM1;

        push!(allmats, DM1);
        push!(allmats2, DM2);
        push!(dampmats, DM2*DM2);
        push!(grids, grid1);
    end

    return allmats2, allmats, grids, dampmats

end



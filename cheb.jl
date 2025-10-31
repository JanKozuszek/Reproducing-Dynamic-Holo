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

function chebmat(grid::AbstractVector{BigFloat})
    # Compute a differentiation matrix given a (Chebyshev) grid

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

    return(Dmat);
end;

function cheb(start, stop, length)
    grid = chebgrid(start, stop, length);
    mat = chebmat(grid);

    return([mat, grid]);
end;

function MultiGridChebyshev(start, stop, Ndom, Npts)
    allmats = zeros(T,Ndom,Npts,Npts);
    allmats2 = zeros(T,Ndom,Npts,Npts);
    grid = zeros(T, Ndom * Npts);
    L = (stop - start)/Ndom;
    loc_start = copy(start);
    loc_stop = loc_start + L;

    for ii in 1:Ndom
        allmats[ii,1:Npts,1:Npts], grid[(ii - 1)*Npts + 1: ii*Npts] = cheb(loc_start, loc_stop, Npts);
        allmats2[ii,1:Npts,1:Npts] = allmats[ii,1:Npts,1:Npts]*allmats[ii,1:Npts,1:Npts];
        loc_start = loc_stop;
        loc_stop = loc_start+L;
    end;

    return allmats2, allmats, grid;
end;
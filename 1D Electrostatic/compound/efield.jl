function efield(rho, dx, Lx, numCells, CALL)
if CALL == 0
    nrho = -(dx^2) .* rho;

    M1 = 1 .* ones(Lx);
    M2 = -2 .* ones(Lx+1);
    M3 = 1 .* ones(Lx);
    M = Tridiagonal(M1, M2, M3);

    phi = M \ nrho;

    ex = zeros(numCells);

    for i = 1:numCells
        # ex[i] = -(phi[i+1] - phi[i-1]) ./ (2*dx);
        if i == 1 || i == numCells
            ex[i] = -(phi[2] - phi[end-1]) ./ (2*dx);
        else
            ex[i] = -(phi[i+1] - phi[i-1]) ./ (2*dx);
        end
    end

    ex[101] = ex[100];
    ex[1] = ex[2];

    return phi, ex;

elseif CALL == 9

    println("STOPSTOPSTOP")

    M1 = ones(Lx);
    M2 = -2 .* ones(Lx+1);
    M3 = ones(Lx);
    M = Tridiagonal(M1, M2, M3);

py""" 
n = $numCells
k = np.array([np.ones(n-1),-2*np.ones(n),np.ones(n-1)])
offset = [-1,0,1]
tridiag = diags(k,offset).toarray()
# dx = $dx
# rho = -(dx**2)*$rho
"""

py"""
phi = -np.linalg.solve(tridiag,-($dx**2)*$rho)
"""

jphi = M \ rho;

pyphi = py"phi";

for i = 1:101
    if pyphi[i] < (jphi[i]-0.1)
        println("pyphi[", i, "] = ", pyphi[i])
        println("jphi[", i, "] = ", jphi[i])
        error("Problem in φ- calculation")
    elseif pyphi[i] > (jphi[i]+0.1)
        println("pyphi[", i, "] = ", pyphi[i])
        println("jphi[", i, "] = ", jphi[i])
        error("Problem in φ+ calculation")
    end
end

phi = pyphi;

py"""
# ex = -np.gradient(phi, $dx, edge_order=2)
ex = np.zeros(numcells)
ex[1] = -(phi[2] - phi[-2])/(2*$dx)
ex[-1] = ex[1]
ex[2:-1] = (phi[:-2].real - phi[2:].real)/(2*$dx)
"""

pyex = py"ex" # = np.gradient($phi, $dx, edge_order=2)"
jex = zeros(numCells);

for i = 1:numCells
    if i == 1 || i == numCells
        jex[i] = -(phi[2] - phi[end-1]) ./ (2*dx);
    else
        jex[i] = -(phi[i+1] - phi[i-1]) ./ (2*dx);
    end
end

for i = 1:101
    if pyex[i] < (jex[i]-0.7)
        println("pyex[", i, "] = ", pyex[i])
        println("jex[", i, "] = ", jex[i])
        error("Problem in ex- calculation")
    elseif pyex[i] > (jex[i]+0.7)
        println("pyex[", i, "] = ", pyex[i])
        println("jex[", i, "] = ", jex[i])
        error("Problem in ex+ calculation")
    end
end

return pyphi, pyex;

elseif CALL == 1
py"""
n = $numCells
k = np.array([np.ones(n-1),-2*np.ones(n),np.ones(n-1)])
offset = [-1,0,1]
tridiag = diags(k,offset).toarray()

phi = np.linalg.solve(tridiag,-($dx**2)*$rho)
ex = -np.gradient(phi, $dx, edge_order=2)
ex = np.zeros(n)
ex[0] = (phi[-1].real - phi[1].real)/(2*$dx)
ex[-1] = ex[0]
# ex[1:-1] = (phi[:-2].real - phi[2:].real)/(2*$dx)
"""

    ex = py"ex";
    phi = py"phi";

    return phi, ex;
end

end
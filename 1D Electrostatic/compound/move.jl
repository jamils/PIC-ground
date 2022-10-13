function move(x, vx, n, dt, q, numCells, dx, Lx, CALL)
    if CALL == 0

        x = x .+ vx .* dt;
        x = mod.(x,Lx);

        rho = zeros(numCells);

        for i = 1:n
            xl = floor(x[i]);
            xu = ceil(x[i]);

            ui = Int64(xu + 1);
            li = Int64(xl + 1);

            rho[li] = rho[li] + ((xu - x[i]) * (q));
            rho[ui] = rho[ui] + ((x[i] - xl) * (q));
        end

        # rho[1] = rho[1] + rho[end];
        # rho[end] = rho[1];

        rho[end] = rho[1] + rho[end];
        rho[1] = rho[end];

        return x, rho;

    elseif CALL == 1
        x = x .+ vx .* dt;
        x = mod.(x,Lx);

py"""
count = 0
rho = np.zeros($numCells)
x = $x;
for j in x:
    fj = np.floor(j)
    cj = np.ceil(j)
    rho[int(fj)] += $q * (cj - j)/$dx
    rho[int(cj)] += $q * (j - fj)/$dx
    count += 1
rho[-1] = rho[0]+rho[-1]
rho[0] = rho[-1]
"""

        rho = py"rho";

        return x, rho;
    end

end
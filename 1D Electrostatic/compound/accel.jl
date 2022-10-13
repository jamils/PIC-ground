function accel(vx, x, n, qm, dt, dx, ex, CALL)
    if CALL == 0
        E = zeros(n);

        for i = 1:n
            xu = ceil(x[i]);
            xl = floor(x[i]);

            ui = Int64(xu + 1); # + 1
            li = Int64(xl + 1); # + 1

            E[i] = ex[li]*((xu - x[i])/dx) + ex[ui]*((x[i] - xl)/dx);
        end
        # println(lineplot(E))

        vx += (qm * dt) .* E;

        return vx;

    elseif CALL == 1
py"""
vxn = np.ones(len($vx))
dx = $dx

count = 0
for j in $x:
    if j>$Lx or j<0:
        j = j %Lx
        $x[count] = j
        E = 0
        fj = np.floor(j)
        cj = np.ceil(j)
        E = ex[int(fj)] * ((cj - j)/dx) + ex[int(cj)] * ((j - fj)/dx)
        vxn[count] = $vx[count] + ($qm * E * $dt)
        count += 1
"""
        # vxn[count] = $vx[count-1] + ($qm * E * $dt)/(2*$me) + ($qm * E * $dt)/(2*$me)

        vx = py"vxn";
        

        return vx
    end

end
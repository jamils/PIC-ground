println("Setting initial constants")

# Initial Conditions

numCells = 101;
Lx = 100;
numSteps = 5000;
dt = 0.05;
dx = Lx/(numCells-1);

# Load particle postions and velocities

ne = 100;
nb = 2000;
v0e = 0;
v0b = 16;
wpe = 1.0;
wpb = 0.032;
qme = -0.001;
qmb = -1.0;

println("Calculating initial conditions")

qe = wpe^2/qme/(ne/Lx);
qb = wpb^2/qmb/(nb/Lx);
me = qe/qme;
mb = qb/qmb;
dpe = 0.5 * Lx/ne;
dpb = 0.5 * Lx/nb;

xe = range(dpe, Lx-dpe, length=ne);
xb = range(dpb, Lx-dpb, length=nb);

xb = xb + 0.001 .* sin.(2 .* pi .* xb ./ Lx);

vxe = v0e .* ones(size(xe));
vxb = v0b .* ones(size(xb));

# Initialize charge density, potential, and e-field

println("Initializing rhoe")

# rhoe = zeros(numCells);
# je = broadcast(floor,xe);
# dle = (xe - je) ./ dx;
# due = (ones(size(dle)) - dle) ./ dx;
# dle = qe .* dle;
# due = qe .* due;

# for i = 1:ne
#     rhoe[Int64(je[i] + 1)] = rhoe[Int64(je[i] + 1)] + due[i];
#     rhoe[Int64(je[i] + 2)] = rhoe[Int64(je[i] + 2)] + dle[i];
# end

# rhoe[end] = rhoe[end] + rhoe[1];
# rhoe[1] = rhoe[end];



println("Initializing rhob")

xe, rhoe = move(xe, vxe, ne, dt, qe, numCells, dx, Lx, 0);
xb, rhob = move(xb, vxb, nb, dt, qb, numCells, dx, Lx, 0);

# rhob = zeros(numCells);
# jb = broadcast(floor,xb);
# dlb = (xb - jb) ./ dx;
# dub = (ones(size(dlb)) - dlb) ./ dx;
# dlb = qb .* dlb;
# dub = qb .* dub;

# for i = 1:nb
#     rhob[Int64(jb[i] + 1)] = rhob[Int64(jb[i] + 1)] + dlb[i];
#     rhob[Int64(jb[i] + 2)] = rhob[Int64(jb[i] + 2)] + dub[i];
# end

# rhob[end] = rhob[end] + rhob[1];
# rhob[1] = rhob[end];

println("Initializing rhoi")

rhoi = ((-qe * ne - qb * nb) / Lx) .* ones(numCells);
println("Calculating rho")
rho = rhoe + rhob + rhoi;

scatterplot(xe)
scatterplot(xb)
scatterplot(rho)
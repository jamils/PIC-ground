using ForwardDiff
using LinearAlgebra

function accel(vx, x, qm, dt, dx, ex, numcells, Lx)
    count = 0;
    for j in x
        if j>Lx || j<0
            j = mod(j,Lx);
            x[count] = j;
        end
        E = (ex[floor(Int,j)] * (ceil(Int,j)-j)/dx) + ex[ceil(Int,j)]*(j-floor(j));
        vx[count] = vx[count-1] + (qm * E * dt)/(2*me) + (qm * E * dt)/(2*me);
        count += 1;
    end

    return vx;
end

function move(x, vx, n, dt, q, numcells, dx, Lx)
    x = x + (vx * dt);
    count = 0;
    rho = zeros(numcells);
    for j in x
        if j>Lx || j<0
            j = mod(j,Lx);
            x[count] = j;
        end
        rho[floor(Int,j)] += (q/Lx) * ((ceil(Int,j)-j)/dx);
        rho[ceil(Int,j)] += (q/Lx) * ((j - floor(Int,j))/dx);
        count += 1
    end

    return x,rho
end

function efield(rho, dx, numcells, Lx, tridiag)
    dx = Lx/(numcells-1);
    phi = tridiag/(-(dx^2)*rho);
    E = ForwardDiff.gradient(phi,dx);
    E[0] = (phi[-1] - phi[1])/(2*dx)
    E[-1] = E[0]

    return phi,E 
end

numcells = 101;
Lx = 100;
numsteps = 5000;
dt = 0.005;
dx = Lx/(numcells-1);

klow = ones(numcells-1);
kmid = -2 * ones(numcells);
khigh = ones(numcells-1);
tridiag = LinearAlgebra.Tridiagonal(klow, kmid, khigh)

# load particle positions and velocities
ne  = 2000;               # number of background electrons %%%%%%%%%%% Should be 100
nb  = 2000;              # number of beam electrons
v0e = -16;                 # background electron drift velocity
v0b = 16;                # beam electron drift velocity
wpe = 0.032;              # background electron plasma frequency
wpb = 0.032;             # beam electron plasma frequency
qme = -1.000;            # background electron charge to mass
                         #(artificially small for linear background motion)
qmb = -1.000;            # beam electron charge to mass
qe  = wpe^2/qme/(ne/Lx); # background electron charge (signed)
qb  = wpb^2/qmb/(nb/Lx); # beam electron charge       (signed)
me  = qe/qme;            # background electron mass
mb  = qb/qmb;            # beam electron mass
dpe = 0.5*Lx/ne;
dpb = 0.5*Lx/nb;

xe  = range(dpe, stop=Lx-dpe, length=ne);       # background electron positions
xb  = range(dpb, stop=Lx-dpb, length=nb);       # beam electron positions

#for i in 1:2000
#    toms = xb[i] + sin(2 * pi * xb[i]/Lx)
#    xb[i] = toms
#end

#xb  = xb + 0.0001*sin(2*3.14159*xb/Lx);  # beam position perturbation
                                       # (which grows Into a plasma wave) 

vxe = v0e*ones(ne);          # background electron velocity
vxb = v0b*ones(nb);          # beam electron velocity

#initialize charge density, potential, and electric field
rhoe = zeros(numcells);
je   = floor.(Int, xe);
dle  = (xe-je)/dx;
due  = (ones(ne) - dle)/dx;
dle  = qe*dle;
due  = qe*due;

rhoe[Int(je-1)+1] = rhoe[Int(je-1)+1] + due;
rhoe[Int(je-1)+2] = rhoe[Int(je-1)+2] + dle;

rhoe[Int(numcells-1)] = rhoe[Int(numcells-1)] + rhoe[0];
rhoe[0] = rhoe[Int(numcells-1)];

rhob = zeros((Int(numcells)));
jb   = floor.(Int,xb);
dlb  = (xb-jb)/dx;
dub  = (ones(nb)-dlb)/dx;
dlb  = qb*dlb;
dub  = qb*dub;

for i in 0:Int(nb):
    rhob[Int(jb[i]-1)+1] = rhob[Int(jb[i]-1)+1] + dlb[i];
    rhob[Int(jb[i]-1)+2] = rhob[Int(jb[i]-1)+2]+ dub[i];
end

rhob[Int(numcells-1)] = rhob[Int(numcells-1)] + rhob[0];
rhob[0] = rhob[Int(numcells-1)];

# assume ions are stationary to collect the ion charge density
rhoi = ((-qe*ne-qb*nb)/Lx)*ones((Int(numcells)));
# total charge density from electrons and ions
rho = rhoe + rhob + rhoi;

phi,ex = efield(rho, dx, numcells, Lx, tridiag);
energy = zeros(numsteps)

for i in 0:numsteps
    vxe = accel(vxe, xe, qme, dt, dx, ex, numcells, Lx)
    vxb = accel(vxb, xb, qmb, dt, dx, ex, numcells, Lx)

    xe,rhoe = move(xe, vxe, ne, dt, qe, numcells, dx, Lx)
    xb,rhob = move(xb, vxb, nb, dt, qb, numcells, dx, Lx)

    rhoe[0] = rhoe[-1] + rhoe[0]
    rhoe[-1] = rhoe[0]

    rhob[0] = rhob[-1] + rhob[0]
    rhob[-1] = rhob[0]

    rho = rhoe + rhob + rhoi

    phi,ex = efield(rho, dx, numcells, Lx, tridiag)

    we = .5 * sum(ex^2)
    kep = .5 * sum(.5*me*(vxe^2))
    keb = .5 * sum(.5*mb*(vxb^2))
    energy[i] = we + kep + keb

    
end
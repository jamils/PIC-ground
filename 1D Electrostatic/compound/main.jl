println("Starting...")

using Pkg
# using PyPlot
using Plots
# using PyCall
using LinearAlgebra
using UnicodePlots
using ProgressMeter
# using Markdown
# using CuArrays

#= py"""
import numpy as np
import math
from scipy.interpolate import interp1d
import scipy.linalg
from scipy.sparse import diags
from numpy.linalg import inv
""" =#

#########
# FLAGS 

GIF = 0     # If you want gifs outputted at the end of the simulation of the evolution of each value. This will significantly increase the simulation time. 1 to enable GIF output, 0 to disable

BEAM = 0    # 0 to set initial conditions for the Weak Beam instability, 1 to set for Two-Stream instability

SECT = 4    # Used to define which python function you want outputted in the GIFs. This follows as the table describes below:

#= 
1 - ACCEL
2 - MOVE
3 - EFIELD
4 - No Tests 
=#

#########

include("efield.jl")
include("move.jl")
include("accel.jl")
include("init.jl")

# ------------------------------------------------------------------------
# ----------------------------- Main Program -----------------------------
# ------------------------------------------------------------------------

# Creating initial arrays to reset in loop
xeinit = xe;
xbinit = xb;
vxeinit = vxe;
vxbinit = vxb;
rhoeinit = rhoe;
rhobinit = rhob;
rhoinit = rho;

# Creating history for arrays
xeh = zeros(100, 4, numSteps);
xbh = zeros(2000, 4, numSteps);
vxeh = zeros(100, 4, numSteps);
vxbh = zeros(2000, 4, numSteps);
rhoeh = zeros(101, 4, numSteps);
rhobh = zeros(101, 4, numSteps);
rhoh = zeros(101, 4, numSteps);
phih = zeros(101, 4, numSteps);
exh = zeros(101, 4, numSteps);

count = 0;

for i = 4:4 # 1:4

    println(i,": ---------------")

    xe = xeinit;
    xb = xbinit;
    vxe = vxeinit;
    vxb = vxbinit;
    rhoe = rhoeinit;
    rhob = rhobinit;
    rho = rhoinit;

    global count = 3;

    if count == 0
        println("Intial electric field computation")
        phi, ex = efield(rho, dx, Lx, numCells, 1);
    else
        println("Intial electric field computation")
        phi, ex = efield(rho, dx, Lx, numCells, 0);
    end
    global phiinit = phi;
    global exinit = ex;
    x = 0:100;

    scatterplot(phi)
    scatterplot(ex)

    # Compute Energies
    global et = zeros(numSteps);
    global ef = zeros(numSteps);
    global keb = zeros(numSteps);
    global kee = zeros(numSteps);

    acceltime = zeros(numSteps);
    movetime = zeros(numSteps);
    efieldtime = zeros(numSteps);

    # ------------------------------------------------------------------------

# md"""
# ## Testing ACCEL subroutine

# This step runs all of the standard Julia functions, with an exception for the ACCEL subroutine - this is ran using the Python code. The general structure of this subroutine is:

# * $$\vec{v} = \vec{v}^{n-\frac{1}{2}} + \frac{q \vec{E}}{m} \frac{\Delta t}{2}$$
# * $$\vec{v}_1 = \vec{v}^- + \vec{v} \times \vec{t}$$
# * $$\vec{v}^+ = \vec{v}^- + \vec{v}_1 \times \vec{s}$$
# * $$\vec{v}^{n + \half} = \vec{v}^+ + \frac{q \vec{E}}{m} \frac{\Delta t}{2}$$

# """

    # ------------------------------------------------------------------------

    println("Main Loop starting...")
    prog = Progress(numSteps, 1)
    for j = 1:numSteps
        if count == 2
            # println("---Testing accel function---")
            vxe = accel(vxe, xe, ne, qme, dt, dx, ex, 1);
            vxb = accel(vxb, xb, nb, qmb, dt, dx, ex, 1);
            # acceltime[j] = @elapsed accel(vxb, xe, ne, qme, dt, dx, ex);

            xe, rhoe = move(xe, vxe, ne, dt, qe, numCells, dx, Lx, 0);
            xb, rhob = move(xb, vxb, nb, dt, qb, numCells, dx, Lx, 0);
            # movetime[j] = @elapsed move(xe, vxb, ne, dt, qe, numCells, dx, Lx);

            rho = rhoe + rhob + rhoi;

            phi, ex = efield(rho, dx, Lx, numCells, 0);
            # efieldtime[j] = @elapsed efield(rho, dx, Lx, numCells);

        elseif count == 1
            # println("---Testing move function---")
            vxe = accel(vxe, xe, ne, qme, dt, dx, ex, 0);
            vxb = accel(vxb, xb, nb, qmb, dt, dx, ex, 0);
            # acceltime[j] = @elapsed accel(vxb, xe, ne, qme, dt, dx, ex);

            xe, rhoe = move(xe, vxe, ne, dt, qe, numCells, dx, Lx, 1);
            xb, rhob = move(xb, vxb, nb, dt, qb, numCells, dx, Lx, 1);
            # movetime[j] = @elapsed move(xe, vxb, ne, dt, qe, numCells, dx, Lx);

            rho = rhoe + rhob + rhoi;

            phi, ex = efield(rho, dx, Lx, numCells, 0);
            # efieldtime[j] = @elapsed efield(rho, dx, Lx, numCells);

        elseif count == 0
            # println("---Testing efield function---")
            vxe = accel(vxe, xe, ne, qme, dt, dx, ex, 0);
            vxb = accel(vxb, xb, nb, qmb, dt, dx, ex, 0);
            # acceltime[j] = @elapsed accel(vxb, xe, ne, qme, dt, dx, ex);

            xe, rhoe = move(xe, vxe, ne, dt, qe, numCells, dx, Lx, 0);
            xb, rhob = move(xb, vxb, nb, dt, qb, numCells, dx, Lx, 0);
            # movetime[j] = @elapsed move(xe, vxb, ne, dt, qe, numCells, dx, Lx);

            rho = rhoe + rhob + rhoi;

            phi, ex = efield(rho, dx, Lx, numCells, 1);
            # efieldtime[j] = @elapsed efield(rho, dx, Lx, numCells);

        elseif count == 3
            # println("---Correct solution---")
            vxe = accel(vxe, xe, ne, qme, dt, dx, ex, 0);
            vxb = accel(vxb, xb, nb, qmb, dt, dx, ex, 0);
            # acceltime[j] = @elapsed accel(vxb, xe, ne, qme, dt, dx, ex);

            xe, rhoe = move(xe, vxe, ne, dt, qe, numCells, dx, Lx, 0);
            xb, rhob = move(xb, vxb, nb, dt, qb, numCells, dx, Lx, 0);
            # movetime[j] = @elapsed move(xe, vxb, ne, dt, qe, numCells, dx, Lx);

            rho = rhoe + rhob + rhoi;

            phi, ex = efield(rho, dx, Lx, numCells, 0);
            # efieldtime[j] = @elapsed efield(rho, dx, Lx, numCells);
        end

        ef[j] = 0.5 * sum(ex.^2);
        kee[j] = sum(0.5 * me * (vxe .^2));
        keb[j] = sum(0.5 * mb * (vxb .^2));
        et[j] = ef[j] + keb[j] + kee[j];
        # return ef, kee, keb, et;

        # Add to histories

        xeh[:, i, j] = xe;
        xbh[:, count+1, j] = xb;
        vxeh[:, count+1, j] = vxe;
        vxbh[:, count+1, j] = vxb;
        rhoeh[:, count+1, j] = rhoe;
        rhobh[:, count+1, j] = rhob;
        rhoh[:, count+1, j] = rho;
        phih[:, count+1, j] = phi;
        exh[:, count+1, j] = ex;

        # Advance progress bar

        next!(prog)
    end

    println("Change in total energy: ", et[end]-et[1])

    count = count + 1;

    # xeh[:, count, :] = xe;
    # xbh[:, count, :] = xb;
    # vxeh[:, count, :] = vxe;
    # vxbh[:, count, :] = vxb;
    # rhoeh[:, count, :] = rhoe;
    # rhobh[:, count, :] = rhob;
    # rhoh[:, count, :] = rho;
    # phih[:, count, :] = phi;
    # exh[:, count, :] = ex;
    # efh[:, count, :] = ef;
    # kebh[:, count, :] = keb;
    # keeh[:, count, :] = kee;

    println(" ")

end

println("Sumlations completed successfully")


###

# Use @gif and markdown at bottom with flags at top, not in main loop

# GIF either on or off, and plots all or none. Need to change history arrays to sort arrays and create simulation history arrays for each variable. 

# Need to rewrite python so that I can analyze each step, and figure out how to have Julia compare values and provide rolling histograms with discrepency bars for each value for each variable for each timestep.
    # Look at each function, find the one with significant difference in values, and then create histograms for that particular fuction with md for algorithm

# Write all of this up, with md contextual information

if GIF == 1
    println("Creating GIF output...")
    progvid = Progress(numSteps, 1)
    anim = @animate for i = 1:numSteps
        scatter(1:2000, vxbh[:, SECT, i])
        # scatter(1:100, vxeh[:, SECT, i])
        # plot([rhoh[:, SECT, i], phih[:, SECT, i]])

        next!(progvid)
    end every 10

    gif(anim, fps = 8) 
else 
    println("No GIF outputted")
end

println("∎")
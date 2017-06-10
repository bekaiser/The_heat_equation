# Bathtub advection / diffusion / latent heat flux problem
# Bryan Kaiser
# 9/25/16

# NOTE: make the directory "/T" for output plots.

using DataArrays
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport pylab as py
using HDF5

# import functions
include("./functions.jl")


# =============================================================================
# physical parameters

# fluid properties 
c = 4.0e6; # J/(m^3*K), equivalent to rho*cp (product of specific heat and density)
h = 3.0e3; # W/(m^2*K), heat flux coefficient
kappa = 5.0e3; # W/(m*K), diffusion coefficient
d = 0.5; # m, tub depth
L = 2.0; # m, tub length
TC = 20.0; # C, (K-273.15), air temperature (cold)
TH = 30.0; # C, (K-273.15), inflow water (hot)

U = 0.004; # m/s, inflowing/outflowing velocity magnitude
K = kappa/c; # m^2/s, diffusivity
H = h/(c*d); # 1/s

# grid
N = 101; dx = L/float(N-1); x = collect(0.0:dx:L);

# analytical solution
Tsoln = bathtub_analytical(c,d,h,kappa,L,TC,TH,U,x);

# initial condition
T0 = x.*0.0+25.0; 
Tn = T0; # T[n], n = time step
Tnm1 = T0; # T[n-1], n-1 = previous time step
Tsave = T0.*0.0;

# time steppin' (run to 2500 s)
dt = 0.1; tn = 0.0; 
# U=0.01: N=21 dt=3.0 | | N=101, dt=0.1 | N=101, dt=0.02
# U=0.001:  N=21 dt=5.0 |N=51 dt=1.0 | N=101, dt=0.2 | N=201, dt=0.05 
dt_2 = dt/2.0;
dt_6 = dt/6.0;
n0 = 0; # starting time step (integer)
Nt = Int32(3000.0/dt); # maximum time step for 2500s

# stability characteristics
Co = U*dt/dx; Di = K*dt/(dx^2.0);
println("\nCourant number = $Co, Diffusion number = $Di\n");

# convergence criteria
dT_sum = 1.0; # initial value
dT_sum_crit = 4E-10; # stop the simulation when the previous time step is only 
# different from the current time step by this amount.

# output information
plotflag = 1; # enter 1 for .png output plots or 0 for no output plots
nplot_start = 0; # number of time steps after n0 before starting output plots
nplot_interval = 100; # number of time steps between output plots
n10 = 0; # counter for time steps
nplot = 0; # counter for output plots
writename = "./sim_N101_U004_w016.h5";


# =============================================================================
# finite difference operator (Pade 4th-order)

N = Int32(N); 

# 1st derivative: 4th-order Pade finite difference operator, dTdx = A1\b1
A = diagm(ones(N).*(4.0*dx/3.0),0)+diagm(ones(N-1).*(dx/3.0),1)+diagm(ones(N-1).*(dx/3.0),-1);
A[1,1:2] = [1.0 0.0]; A[N,N-1:N] = [0.0 1.0]; # for dT/dx, d^2T/dx^2 boundaries
b1 = zeros(N,1); # Pade RHS vector, 1st derivative  
b2 = zeros(N,1); # Pade RHS vector, 2nd derivative  


# =============================================================================
# simulation

n = 0;
tic()
for n = n0:Nt+n0
	
	Tnm1 = Tn; # temperature profile at previous step (for convergence check)
	if n >= 500
        Tsave = Tsave+Tn;
	end

	# terminal output
	if (n-n0)/10 == n10 # show time step every 10 steps
	println("\ntime: $(tn) s"); # time elapsed
	println("$(toq()/10) s per time step"); # wall time / time step
	println("dT_sum = $dT_sum"); # convergence criteria
	tic()
	n10 = n10+1;
	end 

	# output plots for a movie
	if plotflag == 1 
	if (n-n0) >= nplot_start # start output
	if (n-nplot_start-n0)/nplot_interval == nplot # plot interval
	fig1 = py.figure() # psi figsize=(plotsize)
	CP1 = py.plot(x,T0,"r",label="T0");
	CP2 = py.plot(x,Tn,"b",label="T"); 
        py.axis([0.0, 2.0, TC, TH])
	py.xlabel("x"); py.ylabel("T"); py.legend()
	if n <10 # filename
	plotname = "./T/000000$(convert(Int32,n)).png"; elseif n <100
	plotname = "./T/00000$(convert(Int32,n)).png"; elseif n <1000
	plotname = "./T/0000$(convert(Int32,n)).png"; elseif n <10000
	plotname = "./T/000$(convert(Int32,n)).png"; elseif n <100000
	plotname = "./T/00$(convert(Int32,n)).png"; elseif n <1000000
	plotname = "./T/0$(convert(Int32,n)).png"; elseif n <10000000
	plotname = "./T/$(convert(Int32,n)).png";
	end
	py.savefig(plotname,format="png"); py.close(fig1);
	nplot = nplot+1
	end # plot interval
	end # start outputing .png plots into /T folder
	end # plots on (plotflag)

	# 4th-order Runge-Kutta time advance 
	k1 = RK4(A,b1,b2,dx,H,Tn,tn,K,TC,TH,U); # 1st RK4 coefficient
	k2 = RK4(A,b1,b2,dx,H,Tn+k1.*dt_2,tn+dt_2,K,TC,TH,U); # 2nd RK4 coefficient
	k3 = RK4(A,b1,b2,dx,H,Tn+k2.*dt_2,tn+dt_2,K,TC,TH,U); # 3rd RK4 coefficient
	k4 = RK4(A,b1,b2,dx,H,Tn+k3.*dt,tn+dt,K,TC,TH,U); # 4th RK4 coefficient
	Tn = Tn+(k1+k2.*2.0+k3.*2.0+k4).*dt_6; # advance T
	tn = tn+dt; # advance time

	# convergence criteria
	dT_sum = abs(sum(Tn)-sum(Tnm1));
	if dT_sum <= dT_sum_crit	
		break
	end

end # time advancement loop
toq()


# =============================================================================
# post processing

mu_T = Tsave/Float32(Nt) # mean temperature
println("mu_T=$mu_T");

#=

# save results to file
n_conv = n; # converged time step
mass_flux = ((U*tn)*d^2.0)*1000.0; # m*m^3 = m^3*1000 = liters
std_T = std(Tn); # standard deviation (minimize this)
mu_T = mean(Tn); # mean temperature
h5open(writename, "w") do file
    	write(file,"Co",Co)  # alternatively, say "@write file A"
	write(file,"Di",Di); write(file,"c",c); write(file,"h",h)  
	write(file,"kappa",kappa); write(file,"d",d); write(file,"Tsoln",Tsoln);
	write(file,"mass_flux",mass_flux); #write(file,"n_conv",n_conv);
	write(file,"tn",tn); write(file,"Nt",Nt); write(file,"N",N);
	write(file,"L",L); write(file,"TC",TC); write(file,"TH",TH); 
	write(file,"U",U); write(file,"K",K); write(file,"H",H); write(file,"Tsave",Tsave);
	write(file,"std_T",std_T); write(file,"mu_T",mu_T); write(file,"Tn",Tn); 
	#write(file,"dT_sum_crit",dT_sum_crit); write(file,"error",error); 
end # filename

# solution comparison plot
fig2 = py.figure() # psi figsize=(plotsize)
CP1 = py.plot(x,Tsoln,"g",label="analytical");
CP2 = py.plot(x,Tn,"--b",label="computed"); #py.axis([0.0, 2.0, TC, TH])
py.xlabel("x"); py.ylabel("T"); py.legend()
plotname = "./solution.png";
py.savefig(plotname,format="png"); py.close(fig2);

# solution error plot
error = abs(Tsoln)-abs(Tn);
fig3 = py.figure();
CP1 = py.semilogy(x,error,"k",label="analytical");
py.xlabel("x"); py.ylabel("error"); py.title("resolution N=$N");
plotname = "./solution_error.png"; 
py.savefig(plotname,format="png"); py.close(fig3);

=#





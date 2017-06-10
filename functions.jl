# Functions for advection / diffusion / latent heat flux bathtub problem
# Bryan Kaiser
# 9/25/16


# =============================================================================
# function declaration 

function RK4(A,b1,b2,dx,H,T,t,K,TC,TH,u0)
    # fourth-order Runge-Kutta coefficients
    
    #U = u0/2.0; # time independent U 
    U = u0/2.0*(sin(4.0*pi/800.0*t)+1.0); # time dependent U 
    
    # 1st derivative (4th-order Pade)
    b1[2:N-1] = T[3:N]-T[1:N-2];
    b1[1] = U*(T[1]-TH)/K; b1[N] = 0.0; # boundary conditions
    dTdx = A\b1;
 
    # 2nd derivative (4th-order Pade)
    b2[2:N-1] = dTdx[3:N]-dTdx[1:N-2];
    #b2[1] = U*dTdx[1]/K; b2[N] = 0.0; # boundary conditions
    b2[1] = (-(25.0/12.0)*dTdx[1]+4.0*dTdx[2]-3.0*dTdx[3]+(4.0/3.0)*dTdx[4]-(1.0/4.0)*dTdx[5])/dx; 
    b2[N] = ((25.0/12.0)*dTdx[N]-4.0*dTdx[N-1]+3.0*dTdx[N-2]-(4.0/3.0)*dTdx[N-3]+(1.0/4.0)*dTdx[N-4])/dx; 
    
    k = -dTdx.*U+(A\b2).*K+(TC-T).*H; # Runge-Kutta coefficient
    return k # Runge-Kutta coefficient 
end


function lambda(c::Float64,d::Float64,h::Float64,k::Float64,U::Float64)
    # exponents and coefficients for the analytical solution

    l1ml2 = (((U*c)/(2.0*k))^2.0+(4.0*h)/(k*d))^0.5; # lambda_1-lambda_2
    l1 = (U*c)/(2.0*k)+0.5*l1ml2; # lambda_1
    l2 = (U*c)/(2.0*k)-0.5*l1ml2; # lambda_1
    a1 = ((1.0-l1/l2*exp(l1ml2*L)*(U*c-k*l2)/(U*c-k*l1))^(-1.0))*((U*c*(TH-TC))/(U*c-k*l1)); # alpha_1
    a2 = -l1/l2*exp(l1ml2*L)*a1; # alpha_2
    return l1,l2,a1,a2
end


function bathtub_analytical(c::Float64,d::Float64,h::Float64,k::Float64,L::Float64,TC::Float64,TH::Float64,U::Float64,x::Array{Float64,1})
    # analytical solution    

    # lambdas and alphas (analytical solution exponents and coefficients)
    l1,l2,a1,a2 = lambda(c,d,h,k,U);   
    # temperature analytical solution profile
    Tsoln = exp(x.*l1).*a1+exp(x.*l2).*a2+TC; 
    
    return Tsoln
end   



import numpy
from scipy.optimize import minimize

def bem(CLCD, pitch, CHORD, BETA, BLADE, v, rpm, rho, blades):
    # pitch is the pitching angle of the blade in degrees
    # CHORD is the vector of station-wise chords along the blade
    # BETA is the vector of station-wise blade angles due to twist, in degrees
    # BLADE is the vector of radial positions of each blade element
    # v is the entry velocity to the disc, in m/s
    # rpm is the rotation per minute of the prop_test_2
    # rho is atmospheric density in metric units
    # blades is the number of prop blades
    
    element = BLADE[-1]-BLADE[-2]; 
    omega = (rpm/60)*2* numpy.pi;                   # calculate element size and angular velocity

    K1 = 0.5*rho*blades*CHORD;                        # calculate equation constants
    K2 = 0.5*rho*blades*CHORD*BLADE;                # calculate equation constants
    K3 = 4.0*numpy.pi*rho*v**2; 
    K4 = K3*omega/v;                                       # calculate equation constants
    THETA = (pitch+BETA)*numpy.pi/180;                     # calculate theta along the blade
    
    A0 = 0.1*numpy.ones_like(BLADE); 
    B0 = 0.01*numpy.ones_like(BLADE);                    # initialize A0, B0
    
    # calculate cost function and gradient at each element along the blade
    def bemfuncmin(ABs):
      A0 = ABs[:len(BLADE)]; 
      B0 = ABs[len(BLADE):];               # unpack the A0 and B0 vectors
      #A0 = A0'; B0 = B0';                                                       # put into row vector form

      DtDr, DqDr = bemsolve(A0, B0);                                             # solve system of equations using estimated A and B values

      TEM1=DtDr/(K3*BLADE*(1+A0)); 
      TEM2=DqDr/(K4*BLADE**3*(1+A0));          # calculated A and B values from system of equations

      ERRA = 0.5*(A0-TEM1)**2; 
      ERRB = 0.5*(B0-TEM2)**2;                         # error function

      dERRA_dA0 = (A0-TEM1); 
      dERRB_dB0 = (B0-TEM2);                             # derivative of the squared error cost functions

      J = numpy.sum(numpy.concatenate((ERRA, ERRB))); 
      grad = numpy.concatenate((dERRA_dA0, dERRB_dB0));         # convert back to column vectors
    
      return J, grad

    
    # solve lift and torque vectors for given A and B vectors
    def bemsolve(A, B):
      V0=v*(1+A); 
      V2=omega*BLADE*(1-B); 
      VLOCSQU=V0**2+V2**2;                   # calculate velocities
      PHI=numpy.atan2(V0,V2);                                        # calculate angles

      CL = CLCD[0]
      CD = CLCD[1]  # get lift and drag coefficients from our lift function

      DtDr=K1*VLOCSQU*(CL*numpy.cos(PHI)-CD*numpy.sin(PHI));                            # calculate thrust per element
      DqDr=K2*VLOCSQU*(CD*numpy.cos(PHI)+CL*numpy.sin(PHI));                     # calculate torque per element
    
      return DtDr, DqDr
    

    costfunction = lambda P: bemfuncmin(P);                                           # establish cost function
    options = {'MaxIterations': 100, 'SpecifyObjectiveGradient': True, "Display": "off"};     # set options
    ABs = numpy.concatenate((A0, B0))                                                        # pack A0, B0 into single column vector

    result = minimize(costfunction, ABs, method='BFGS', jac=True, options=options);                           # solving using unconstrained fmin  
    ABs = result.x

    A0 = ABs[:len(BLADE)]; 
    B0 = ABs[len(BLADE):];               # unpack the A0 and B0 vectors
    #A0 = A0'; B0 = B0';                                                         # put into row vector form

    DT, DQ = bemsolve(A0, B0);                                                 # solve system using converged A0 and B0 values
    thrust=numpy.sum(DT*element); 
    torque=numpy.sum(DQ*element);                                                     # calculate thrust, torque 
    power = torque*omega;                                                       # calculate power

    return thrust, torque, power
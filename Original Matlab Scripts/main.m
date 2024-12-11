function [eff, pwr] = main(m, p, t, hyperparam, param)

    fluidDensity = hyperparam(1);
    fanRPM = hyperparam(2);
    fluidViscosity = hyperparam(3);
    initialPressure = hyperparam(4);
    initialVelocity = hyperparam(5);

    %fan paramters
    chordLength = param(1);
    chordAngle = param(2);
    id = param(3);
    od = param(4);
    bladeNum = param(5);
 
    
%%%%% PHYSICS %%%%%%%

    coords = nacaCoordinates_new(m, p, t);

    X = coords(:, 1);
    Y = coords(:, 2);
    RE=6e6;
    MACH=0;
    alpha=0;

    [CL, CD, converge]=xfoil(X,Y,alpha,RE,MACH);

    CLCD = [CL, CD];

    if (converge == 0)
        eff = -1;
        pwr = NaN;
        return;
    end

    %%% INPUTS FOR BEMT %%%

    chord = 0.012;
    stations = 10;
    element = id/2/stations;
    BLADE = element/2:element:(id/2-element/2);
    BETA = 15*ones(size(BLADE));
    CHORD = chord*ones(size(BLADE));

    avionics = 10;
    dt = 0.1;
    mass = 1.2;
    gravity = 9.81;
    Aeff = 0.11;
    dia = id;

    n = fanRPM/60;

    v = initialVelocity;

    [thrust,torque,power] = bem(CLCD, chordAngle, CHORD, BETA, BLADE, v, fanRPM, fluidDensity, bladeNum);

    thr = 4*thrust;
    DRAG = Aeff*CD*0.5*fluidDensity*v^2;
    dvdt = (thr-DRAG-mass*gravity)/mass;
    v = v+dvdt*dt;

    pwr = 4*power+mass*gravity*v+DRAG*v+avionics;

    J = v/n/dia;
    kt=thrust/(fluidDensity*n*n*dia*dia*dia*dia);
    kq=torque/(fluidDensity*n*n*dia*dia*dia*dia*dia);
    eff = (J/2/pi)*kt/kq;


    %eff = main(bladeParam(1),bladeParam(2), bladeParam(3));
%%%%%%% END PHYSICS %%%%%%%
    %disp(eff);














% coords = nacaCoordinates_new(m,p,t);
% X = coords(:,1);
% Y = coords(:,2);
% RE=6e6;
% MACH=0;
% alpha=0;
% tic
% [p, converge]=xfoil(X,Y,alpha,RE,MACH);
% toc
% if converge == 1
%     CL=p.cl;
%     CD=p.cd;
%     CLCD = [CL, CD];
% else
%     eff = -1;
%     return;
% end
% rad = 0.085;
% chord = 0.012;
% stations = 10;
% element = rad/stations;
% BLADE = element/2:element:(rad-element/2);
% BETA = 15*ones(size(BLADE));
% CHORD = chord*ones(size(BLADE));
% pitch = 0;
% rpm_max = 15000;
% rho = 1.225;
% rho_isa = 1.225;
% blades = 4;
% max_discharge = 0.85;
% avionics = 10;
% dt = 0.1;
% mass = 1.2;
% gravity = 9.81;
% max_alt = 500;
% Aeff = 0.11;
% dia = 2*rad;
% disk_area = pi*rad^2;                                                                         % disk area of a single rotor
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% no_runs = 5;
% 
% for i=0:(no_runs-1)
%     t = 0;
%     v = 0.0001;
%     x = 0.0001;
%     TIME = [];
%     XS = [];
%     VS = [];
%     THRUST = [];
%     POWER = [];
%     EFFICIENCY = [];
%     rpm = rpm_max-i*250;
%     n = rpm/60;
%     while x<max_alt
%         [thrust,torque,power] = bem(CLCD, pitch, CHORD, BETA, BLADE, v, rpm, rho, blades);
%         thr = 4*thrust;
%         DRAG = Aeff*CD*0.5*rho*v^2;
%         dvdt = (thr-DRAG-mass*gravity)/mass;
%         v = v+dvdt*dt;
%         x = x+v*dt;
%         XS = [XS, x];
%         VS = [VS, v];
%         pwr = 4*power+mass*gravity*v+DRAG*v+avionics;
%         THRUST = [THRUST, thr];
%         POWER = [POWER, pwr];
%         TIME = [TIME, t];
%         J = v/n/dia;
%         kt=thrust/(rho*n*n*dia*dia*dia*dia);
%         kq=torque/(rho*n*n*dia*dia*dia*dia*dia);
%         EFFICIENCY=[EFFICIENCY, (J/2/pi)*kt/kq];
%         t = t+dt;
%     end
% end
% 
% eff = mean(EFFICIENCY)
end
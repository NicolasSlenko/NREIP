import numpy
import matplotlib.pyplot as plt

def nacaCoordinates_new(m, p, t):
    n = None;
    recentFoil = None;

    if not(n):
        n = 0;
    
    n = n + 20;
    if (n > 255):
        n = n - 255;

    # Input
    m = 0.095*m;
    p = (0.8*p) + 0.1;
    t = (0.39*t) + 0.01;

    # Number of points
    num_points = 101;
    
    # Generate airfoil coordinates
    coords = generate_naca(m, p, t, num_points);
    
    # Plot airfoil shape
    plt.figure("Naca Foil");

    if not (recentFoil):
        n = 0;
    else:
        del(recentFoil);
    
    plt.figure(1);
    plt.plot(coords[: ,0], coords[:,1], "Color", numpy.random.rand(3));
    #hold on;
    
    recentFoil = plt.plot(coords[: ,0], coords[: ,1], 'r*-');
    
    #for plotting the coordinates
    #axis equal;
    #grid on;
   
    #xlabel('x');
    #ylabel('y');

def generate_naca(m, p, t, num_points):
    
    # Generate x coordinates
    x = numpy.linspace(0, 1, num_points);
    
    # Calculate thickness distribution
    yt = 5 * t * (0.2969 * numpy.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4);
    
    # Calculate camber line
    yc = numpy.zeros(numpy.size(x));
    dyc_dx = numpy.zeros(numpy.size(x));
    for i in numpy.shape(x):
        if (x(i) < p):
            yc(i) = m / (p**2) * (2 * p * x(i) - x(i)**2);
            dyc_dx(i) = (2 * m / p**2) * (p - x(i));
        else:
            yc(i) = m / ((1 - p)**2) * ((1 - 2 * p) + 2 * p * x(i) - x(i)**2);
            dyc_dx(i) = (2 * m / (1 - p)**2) * (p - x(i));
        #end
    #end
    
    # Calculate theta angle
    theta = numpy.atan(dyc_dx);
    
    # Calculate upper and lower surfaces
    xu = x - yt * numpy.sin(theta);
    yu = yc + yt * numpy.cos(theta);
    xl = x + yt * numpy.sin(theta);
    yl = yc - yt * numpy.cos(theta);
    
    # Combine coordinates
    x_coords = (numpy.flip(xu), xl[1, numpy.len-1]);
    y_coords = (numpy.flip(yu), yl[1, numpy.len-1]);
    
     # Ensure even number of rows
    if numpy.mod(numpy.shape(x_coords), 1) != 0: 
        x_coords(numpy.len-1) = [];  # Remove last point to make it even
        y_coords(numpy.len-1) = [];

    # Combine into a single matrix
    coords = numpy.concatenate(x_coords, y_coords);
    return coords;

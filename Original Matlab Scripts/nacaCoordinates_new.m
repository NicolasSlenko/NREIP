function coords = nacaCoordinates_new(m, p, t)
    persistent n;
    persistent recentFoil;
    if (isempty(n))
        n = 0;
    end
    
    n = n + 20;
    if n > 255
        n = n - 255;
    end


    % Input
    m = 0.095*m;
    p = (0.8*p) + 0.1;
    t = (0.39*t) + 0.01;

    % Number of points
    num_points = 101;
    
    % Generate airfoil coordinates
    coords = generate_naca(m, p, t, num_points);
    
    % Plot airfoil shape
    %figure("Naca Foil");

    if (isempty(recentFoil))
        n = 0;
    else
        delete(recentFoil)
    end

    figure(1);
    plot(coords(:,1), coords(:,2), "Color", [rand(1), rand(1), rand(1)]);
    hold on;
    
    recentFoil = plot(coords(:,1), coords(:,2), 'r*-');
    
    axis equal;
    grid on;
   
    xlabel('x');
    ylabel('y');
end

function coords = generate_naca(m, p, t, num_points)
    
    % Generate x coordinates
    x = linspace(0, 1, num_points);
    
    % Calculate thickness distribution
    yt = 5 * t * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x.^2 + 0.2843 * x.^3 - 0.1015 * x.^4);
    
    % Calculate camber line
    yc = zeros(size(x));
    dyc_dx = zeros(size(x));
    for i = 1:length(x)
        if x(i) < p
            yc(i) = m / (p^2) * (2 * p * x(i) - x(i)^2);
            dyc_dx(i) = (2 * m / p^2) * (p - x(i));
        else
            yc(i) = m / ((1 - p)^2) * ((1 - 2 * p) + 2 * p * x(i) - x(i)^2);
            dyc_dx(i) = (2 * m / (1 - p)^2) * (p - x(i));
        end
    end
    
    % Calculate theta angle
    theta = atan(dyc_dx);
    
    % Calculate upper and lower surfaces
    xu = x - yt .* sin(theta);
    yu = yc + yt .* cos(theta);
    xl = x + yt .* sin(theta);
    yl = yc - yt .* cos(theta);
    
    % Combine coordinates
    x_coords = [flip(xu), xl(2:end)];
    y_coords = [flip(yu), yl(2:end)];
    
     % Ensure even number of rows
    if mod(length(x_coords), 2) ~= 0
        x_coords(end) = [];  % Remove last point to make it even
        y_coords(end) = [];
    end

    % Combine into a single matrix
    coords = [x_coords(:), y_coords(:)];
end

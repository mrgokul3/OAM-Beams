function LG = simulateCustomLGBeam(l, p, z, w0, lambda, SP)
    k = (2 * 3.14) / lambda; % wave number
    z_R = (3.14 * (w0^2)) / lambda; % Rayleigh range of the beam

    % Define grid points
    X = linspace(-1.5,1.5, SP);
    Y = linspace(-1.5, 1.5, SP);

    [x, y] = meshgrid(X, Y);

    wz = w0 * sqrt(1 + (z / z_R).^2); % Gaussian beam width
    Rz = z .* (1 + (z_R ./ z).^2); % beam curvature
    phiz = atan(z ./ z_R); % Gouy phase.

    R = sqrt(x.^2 + y.^2); % radial distance from the beam axis
    phi = atan2(y, x); % azimuthal angle

    % The electric field of an LG mode with a topological charge l and radial index p in a cylindrical system coordinates (R, phi, z)
    
    LG = (1 ./ wz) .* (sqrt((2 * factorial(p) / (3.14 * (factorial(abs(l) + p)))))) .* (exp((1i * ((2 * p) + (abs(l) + 1))) * (phiz))) .* (((sqrt(2) * R) ./ wz).^(abs(l))) .* Laguerre(p, abs(l), (2 * R.^2) ./ (wz.^2)) .* exp(-((1i * k * (R.^2)) ./ (2 * Rz))) .* exp(-((R.^2) ./ (wz.^2))) .* exp(1i * l * phi);
    LG = LG ./ max(abs(LG(:)));
    
    end


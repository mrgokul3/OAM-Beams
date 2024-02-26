function HG = hermiteGaussianBeam(nx, ny, lambda, z, SPx, SPy, w0)
    % Define parameters
    k = (2 * 3.14) / lambda; % wave number
    z_R = (3.14 * (w0^2)) / lambda; % Rayleigh range of the beam
    
    % Define grid points
    X = linspace(-5 * w0, 5 * w0, SPx);
    Y = linspace(-5 * w0, 5 * w0, SPy);
    [x, y] = meshgrid(X, Y);
    
    wz = w0 * sqrt(1 + (z / z_R)^2); % Gaussian beam width
    Rz = z * (1 + (z_R / z)^2); % beam curvature
    phiz = atan(z / z_R); % Gouy phase.

    % Hermite Gaussian Beam Simulation

    % In a Cartesian coordinate system (x; y; z), the electric field of a Hermite Gaussian beam can be written as
    HG = (sqrt((2 ^ nx * 2 ^ ny) / (3.14 * factorial(nx) * factorial(ny) * w0^2))) ...
        .* (exp(-1i * k * (x.^2 + y.^2) / (2 * Rz))) ...
        .* (exp(1i * (nx + ny + 1) * phiz)) ...
        .* (exp(-(x.^2 + y.^2) / wz^2)) ...
        .* myHermite(nx, sqrt(2) * x / wz) ...
        .* myHermite(ny, sqrt(2) * y / wz);

    HG = HG ./ max(max(HG));

    function H = myHermite(n, x)
    if n == 0
        H = ones(size(x));
    elseif n == 1
        H = 2 * x;
    else
        H_0 = ones(size(x)); % Hermite polynomial for n = 0
        H_1 = 2 * x; % Hermite polynomial for n = 1
        for k = 2:n
            H = 2 * x .* H_1 - 2 * (k - 1) * H_0;
            H_0 = H_1;
            H_1 = H;
        end
    end
end
end

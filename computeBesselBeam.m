function SP_Bessel = computeBesselBeam(l, z, SP, w0, lambda, xm)
    SPx = SP; % number of sampling points
    SPy = SP; % number of sampling points
    k = (2 * 3.14) / lambda; % wave number
    z_R = (3.14 * (w0 ^ 2)) / lambda; % Rayleigh range of the beam

    % Define grid points
    X = linspace(-5 * w0, 5 * w0, SPx);
    Y = linspace(-5 * w0, 5 * w0, SPy);
    [x, y] = meshgrid(X, Y);

    wz = w0 * sqrt(1 + (z / z_R) ^ 2); % Gaussian beam width
    Rz = z * (1 + (z_R / z) ^ 2); % beam curvature
    phiz = atan(z / z_R); % Gouy phase.

    % Bessel beam simulation
    theta_degrees = 14; % example angle in degrees
    theta_radians = theta_degrees * 3.14 / 180; % convert to radians
    kr = k * sin(theta_radians); % transverse Wavenumber
    kz = k * cos(theta_radians); % transverse Wavenumber
    A0 = 1; % amplitude of the electric field on the optical axis of Bessel beam

    x3 = linspace(-xm, xm, SPx);
    y3 = linspace(-xm, xm, SPy);
    [X3, Y3] = meshgrid(x3, y3);
    R3 = sqrt(X3 .^ 2 + Y3 .^ 2); % radial distance from the beam axis
    phi = atan2(Y3, X3); % azimuthal angle
    J = besselj(l, ((z_R * kr * R3) / (z_R - (1i * z))));
    Phi = exp((1i * l * phi) - (1i * kz * z)) .* exp(((1i * (kr ^ 2) * z * (w0 ^ 2)) - (2 * (k * R3 ^ 2))) / (4 * (z_R - 1i * z)));
    B = A0 * J .* Phi;
    B = B / max(abs(B(:)));

    SP_Bessel = B;
end

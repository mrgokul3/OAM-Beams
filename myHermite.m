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


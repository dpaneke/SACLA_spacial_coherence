tic;

%load('C:\Users\dpane\Documents\Matlab_files\allimportant\initial_values.mat');

xx0 = linspace(-155, 155, 3000);
xx0a = linspace(-1500, 1500, 3000);


Wa = zeros(size(xx0a, 2), size(xx0a, 2));
W = zeros(size(xx0, 2), size(xx0, 2));

for j1 = 1:size(xx0, 2)
    for j2 = 1:size(xx0, 2)
        W(j1, j2) = A^2*exp(-(xx0(j1)^2+xx0(j2)^2)/(4*sig^2))*exp(-(xx0(j1)-xx0(j2))^2/(2*ksi^2));
    end
end

for j1 = 1:size(xx0a, 2)
    for j2 = 1:size(xx0a, 2)
        Wa(j1, j2) = Aa^2*exp(-(xx0(j1)^2+xx0(j2)^2)/(4*siga^2))*exp(-(xx0(j1)-xx0(j2))^2/(2*ksia^2));
    end
end

I = zeros(1, size(xx, 2));
Ia = zeros(1, size(xx, 2));

[X1, X2] = meshgrid(xx0, xx0);

for j = 1:size(xx, 2)
    I(j) = 1/(lam*z0)*trapz(xx0, trapz(xx0, W.*...
    exp(1i*k/(2*z)*((xx(j) - X2).^2 - (xx(j) - X1).^2)), 2));

    percentage(j, size(xx, 2));
end

[X1, X2] = meshgrid(xx0a, xx0a);

for j = 1:size(xx, 2)
    Ia(j) = 1/(lam*z0)*trapz(xx0a, trapz(xx0a, Wa.*...
    exp(1i*k/(2*z0)*((xx(j) - X2).^2 - (xx(j) - X1).^2)), 2));

    percentage(j, size(xx, 2));
end

toc;
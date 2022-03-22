tic;

format long

% block of constants
% distances with um
lam = 1.24e-4; % 0.12 nm
w = 350;       % source half-width
x1 = 2.517e8; % -251.7 metres - source
x2 = 8.5e6;  % 8.45 metres - detector
a = -80;     % (a,b) - graph viewport
b = 370;
slit_down = -55; % slit down position (-21 21) 100 half-width
slit_up =345;    % slit up position
dy = 0.5;          % a step for graph viewport

delta = 0.01;

%define a source
E_g = @(y) exp( -(y./w).^2); 

%define E in the slit plane
E_before_slit = @(y) (1 - 1i)/sqrt(2*lam *x1)*exp(1i*2*pi*x1/lam)*...
integral(@(y1) E_g(y1).*exp(1i*pi*(y1 - y).^2/(lam*x1)), -Inf, Inf);

%define E in the detector plane
E_detector = @(y) -1i/(lam*sqrt(x1*x2))*exp(1i*2*pi*(x1 + x2)/lam)*...
integral2(@(y1, y2) E_g(y1).*exp(1i*pi*((y2 - y1).^2/x1 + (y2 - y).^2/x2)./lam), -Inf, Inf, slit_down, slit_up);

%preparing for plotting
yy = linspace(a,b,(b-a)/dy+1); 

E_g_arr = zeros(size(yy));
E_before_slit_arr = zeros(size(yy));
E_detector_arr = zeros(size(yy));
index = 0;

for y_k = yy
    index = index + 1;
    E_g_arr(index) = E_g(y_k);
    E_before_slit_arr(index) = E_before_slit(y_k);
    E_detector_arr(index) = E_detector(y_k);
end

%plot(exp100yy-195, 0.01.*exp100-0.02, yy, E_detector_arr.*conj(E_detector_arr));
%plot(exp290yy-70, 0.01.*exp290-4e-3, yy, E_detector_arr.*conj(E_detector_arr));
%plot(exp370yy-70, 0.01.*exp370-4e-3, yy, E_detector_arr.*conj(E_detector_arr));

plot(yy, E_detector_arr.*conj(E_detector_arr));

toc;
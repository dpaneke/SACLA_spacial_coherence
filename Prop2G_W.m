function result = Prop2G_W(xx0, xx1, kef, dlt)
    lam = 1.24e-4; % 0.12 nm wavelength 1.24e-4
    z0 = 8.3e6;
    k = 2*pi/lam;
    z = 250e6;
    
    a1 = kef(1); b1 = kef(2); c1 = kef(3);
    a2 = kef(4); b2 = kef(5); c2 = kef(6);
    
    %sig = c1/sqrt(2);
    x1x2 = find_w(kef);
    w = x1x2(2) - x1x2(1); b_eff =(x1x2(1)+x1x2(2))/2;
    
    sig_t = w/(2*sqrt(2*log(2)));
    
    %dlt = 1/sqrt((1/(2*sig_t)^2 + 1/ksi^2));
    
    ksi = 1/sqrt(1/dlt^2 - 1/(2*sig_t)^2);
    zeff1 = k*sig_t*dlt;
    
    q = zeff1/(2*z); % q must be more than 1 or eq
    p = q*(1 + sqrt(1 - 1/q^2));
    
    R = z*(1 + 1/p^2); % D = sqrt(1 + p^2)
   
    I0x = @(x) a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2);
    
    [X1, X2] = meshgrid(xx0, xx0);
    
    %W0 = sqrt(I0x(X1).*I0x(X2)).*exp(-(X2 - X1).^2/(2*ksi^2)).*...
    %    exp(1i*k/(2*R)*((X2-b_eff).^2 - (X1-b_eff).^2));
    
    W0 = I0x((X1+X2)/2).*exp(-(X2 - X1).^2/(2*ksi^2)).*...
        exp(1i*k/(2*R)*((X2-b_eff).^2 - (X1-b_eff).^2));
    
    Iz = zeros(1, size(xx1, 2));
    Fxy = zeros(size(xx0, 2));

    for j = 1:size(xx1, 2)
        Fxy = k/(2*pi*z0)*W0.*exp(1i*k/(2*z0) *...
                    ((X2'-xx1(j)).^2 - (X1'-xx1(j)).^2) );

        Iz(j) = trapz2(xx0, xx0, Fxy);
        percentage(j, size(xx1, 2));
    end
    
    result = Iz;
end
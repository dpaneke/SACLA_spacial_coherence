function result = Prop2G(xx0, xx1, kef, ksi)
    lam = 1.24e-4; % 0.12 nm wavelength 1.24e-4
    z0 = 8.3e6;
    k = 2*pi/lam;
    z = 250e6;
    
    a1 = kef(1); b1 = kef(2); c1 = kef(3);
    a2 = kef(4); b2 = kef(5); c2 = kef(6);
    
    sig = [c1/sqrt(2), c2/sqrt(2)];

    I0 = [a1, a2];

    delt = 1./sqrt(1./(2*sig).^2 + 1/ksi^2);
    zeff1 = k*sig.*delt;

    ze0 = 2*z^2./(zeff1.*(1 + sqrt(1 - (2*z./zeff1).^2)));
    R = z*(1 + (ze0/z).^2);

    bt0 = sqrt(8*pi)*I0.*sig.*delt./(2*sig+delt);
    kp = (2*sig-delt)./(2*sig+delt);

    num_mode = max(1+fix(-2*log(10)./log(kp)));
    %num_mode = 50;
    %disp(num_mode);

    beta1 = zeros(1, num_mode+1);
    beta2 = zeros(1, num_mode+1);

    Em01 = zeros(num_mode+1, size(xx0, 2));
    Em02 = zeros(num_mode+1, size(xx0, 2));

    phi1 = zeros(num_mode+1, size(xx0, 2));
    phi2 = zeros(num_mode+1, size(xx0, 2));

    phi0 = k*z - 2*pi * fix(k*z/(2*pi));

    for ind = 0:num_mode
       beta1(ind+1) = bt0(1)*kp(1)^ind;
       beta2(ind+1) = bt0(2)*kp(2)^ind;

       phi1(ind+1,:) = phi0 - (ind + 1)*atan(1/ze0(1)) + k*(xx0 - b1).^2/(2*R(1));
       phi2(ind+1,:) = phi0 - (ind + 1)*atan(1/ze0(2)) + k*(xx0 - b2).^2/(2*R(2));

       Em01(ind+1,:) = Emode(ind, k, zeff1(1), xx0-b1).*exp(1i*phi1(ind+1,:));
       Em02(ind+1,:) = Emode(ind, k, zeff1(2), xx0-b2).*exp(1i*phi2(ind+1,:));
    end

%     Emz1 = zeros(num_mode+1, size(xx1, 2));
%     Emz2 = zeros(num_mode+1, size(xx1, 2));

%     for ind = 0:num_mode
%        Emz1(ind+1,:) = FresT(xx1, xx0, Em01(ind+1,:), pi/(lam*z0));
%        Emz2(ind+1,:) = FresT(xx1, xx0, Em02(ind+1,:), pi/(lam*z0));
%     end

    Imodes1 = beta1*abs(Em01).^2;
    Imodes2 = beta2*abs(Em02).^2;

    result = [Imodes1; Imodes2];
end


ns = 1.0e-9;          % nanosecond 
mm = 1.0e-3;          % millimeter 
nA = 1.0e-9;          % nano Amperes
mW = 1.0e-3;          % milliWatts
uW = 1.0e-6;          % microWatts

s.h       = 6.626e-34;  % m.^2kg/s, Plank
s.kb      = 1.38e-23;   % Boltzman constant
s.q       = 1.6e-19;    % C, elementary charge
s.c       = 3.0e8;      % m/s, speed of light
s.lambda_ = 1550*ns;    % m, wavelength
s.T0      = 300;        % K, temperature
s.gamma = .1;            % m^-1, medium extinction coeffcient placeholder 

s.P0   = 1*mW;   % W, initial power

s.rho = .1;      % reflectance target
s.D   = 20*mm;   % diameter aperture

s.R0 = 50;        % ohms, photodiode load
s.BW = 1/(1*ns);  % Hz bandwidth
s.Pb = 10*uW;     % background power @ wavelength

s.Lair    = 0.1/100; % scattering absorption losses per per meter
s.Rtarget = 0.1;     % reflectivity of target
s.Tos     = 0.98;    % transmission through optical system

% PD InGaAs
s.Resp = 0.9;                              % A/W,  responsivity
s.Id = 10*nA;                             % dark noise
s.eta = s.h.*s.c/(s.q.*s.lambda_).*s.Resp; % quantum efficiency

% APD InGaAs
s.M   = 10;       % avalanche mult. factor
s.Ids = 10*nA;    % A, surface leakage current
s.Idb = 10*nA;    % A, bulk leakage current
s.F   = 5.5;      % F, excess noise factor (Mcintyre model)

R    = linspace(1,100,1000);


%%
close all
figCreate();
loglog(R,returnPower(s, R)/uW, 'LineWidth',2);
xlabel( "distance R, [m]");
ylabel ("Receiver Optical Power P_r(R), [\mu W]");
title("Optical Power Budget:"+newline+"Target Return Optical Power at the Colocated Detector")
grid('on')
grid('minor')
figSave("OpticaPowerVersusTargetDistance")

%%
figCreate();
semilogx(R, dB(avalachePhotodiodeSNR(s, R)), 'LineWidth',2);
hold("on");
semilogx(R, dB(photodiodeSNR(s, R)), 'LineWidth',2);
xlabel("distance R, [m]");
ylabel ("Photodiode SNR (R), [dB]");
title("Optical Power Budget:"+newline+"Photodiodes SNR vs Target Distance")
legend("APD","PD")
grid('on')
grid('minor')
figSave("PhotodiodesSNRVersusTargetDistance")



%%
function T = transmissionMedium(S, R)
%Beer's law, glug glug
  T = exp(-S.gamma*R);
end
%%
function Pr = returnPower(S, R)
    A  = pi.*(S.D/2).^2;
    Pr = S.P0 .*S.Tos .* transmissionMedium(S, R) .* S.Rtarget./(pi*R.^2)  .* transmissionMedium(S, R) .*A; 
end

%%
function PDSNR = photodiodeSNR(S,R)

    Popt = returnPower(S,R);

    Num = Popt.*S.Resp;

    Den = 2.*S.q.*S.BW.*(Popt.*S.Resp + S.Id) +  4.*S.kb.*S.T0.*S.BW/S.R0;

    PDSNR = Num./Den;
end

%%
function APDSNR = avalachePhotodiodeSNR(S,R)

    Popt = returnPower(S,R);

    Num = Popt.*S.Resp.*S.M;

    Den = 2.*S.q.*S.BW.*(S.Ids + (S.Idb + Popt.*S.Resp).*S.F.*S.M.^2) +  4.*S.kb.*S.T0.*S.BW./S.R0;

    APDSNR =  Num./Den;
end

%%
function v = dB(x)
    v = 10.*log10(x);
end

fprintf("%f %f\n", returnPower(s, 10), photodiodeSNR(s, 10))

%%
function fig = figCreate()
    fig = figure;
    fig.Position = [100 100 1400 800];
end
%%
function figSave(fileLocation)
    saveas(gcf(), fileLocation, 'png')
end    

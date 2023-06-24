%% Example script
%
% This demonstrates how to use the optical-bloch-solver classes to solve
% for the phase shifts involved in Faraday rotation
%
op = opticalSystem('Rb87','D2');                %Sets the optical system to be for Rb87 atoms and the D2 transition
%
% Set laser properties
%
op.laser1.setGaussBeam(100e-6,5e-3)...         %Sets the intensity for a Gaussian beam of 1.6 uW power and 90 um waist
    .setPolarization([1,0,0],'linear')...       %Sets the polarization in the linear basis
    .setStates([2,2],[3,0],3e9);             %Sets the states that the laser addresses. 
%
% Set magnetic field properties. Angles are relative to the direction of
% the laser field (assumed to be in the positive z direction). So th = 0
% and ph = 0 is a magnetic field along the z axis
%
th = 0;
ph = 0;
op.setMagneticField(0.5,...                  %Magnetic field in Gauss
    [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]); %Direction of field
%
% Integration properties. States are ordered in terms of increasing energy,
% so |F=1,mF=1> is 1, |F=2,mF=0> is 2, etc 
%
op.initPop(8) = 1;
op.integrate(1e-7,100e-6);    %Integrate using time step (first argument) up to a given time (second argument) assuming constant fields

%% Plot
%
% Plots the ground state populations
%
figure(1);clf;
op.plotPopulations('ground');
%
% Plots the photon scattering rates per atom
%
figure(2);clf;
op.plotScatteringRates;
grid on;
xlabel('Time [s]');
ylabel('Photon scattering rate [s^{-1}]');
%Display the total number of photons scattered
fprintf(1,'Total Photons = %.2e\n',trapz(op.t,sum(op.getScatteringRates,1)));
%
% Calculates  and plots the absorbance and phase shift for the different
% components of the light field
%
figure(3);clf;
basis = 'linear';
P = op.getPolarization(basis);
if strcmpi(basis,'spherical')
    str = {'-1','0','+1'};
    pol = laser.sphPolBasis*op.laser1.pol;
else
    str = {'x','y','z'};
    pol = op.laser1.pol;
end
for nn = 1:numel(pol)
    if pol(nn) == 0
        P(nn,:) = 0;
    else
        P(nn,:) = P(nn,:)./(op.laser1.field*pol(nn));
    end
end

subplot(2,1,1);
plot(op.t,real(P),'-','linewidth',2);
legend(str);
xlabel('Time [s]');
ylabel('Absorbance');
title('Absorbance');
subplot(2,1,2);
plot(op.t,imag(P),'-','linewidth',2);
legend(str);
xlabel('Time [s]');
ylabel('Reactance');
title('Reactance');
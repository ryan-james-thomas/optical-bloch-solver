%% Example script
%
% This demonstrates how to use the optical-bloch-solver classes to solve
% optical pumping problems
%
op = opticalSystem('Rb87','D2');                    %Sets the optical system to be for Rb87 atoms and the D2 transition
%
% Set laser properties
%
op.laser1.setGaussBeam(10e-3,1e-3)...                    %Sets the intensity in W/m^2. Can also use setGaussBeam(Power [W],waist [m])
    .setPolarization([0,1,0],'spherical')...    %Sets the polarization for [sigma minus, pi, sigma plus] polarizations
    .setStates([1,0],[1,0],00e6);                  %Sets the states that the laser addresses. 
% op.laser2.setGaussBeam(0e-3,10e-3)...                    %Sets the intensity in W/m^2. Can also use setGaussBeam(Power [W],waist [m])
%     .setPolarization([0,1,0],'spherical')...    %Sets the polarization for [sigma minus, pi, sigma plus] polarizations
%     .setStates([2,0],[2,0],00e6);                  %Sets the states that the laser addresses. 
%
% Set magnetic field properties. Angles are relative to the direction of
% the laser field (assumed to be in the positive z direction). So th = 0
% and ph = 0 is a magnetic field along the z axis
%
th = 0*pi/180;
ph = 0;
op.setMagneticField(200e-3,...                  %Magnetic field in Gauss
    [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]); %Direction of field
%
% Integration properties. States are ordered in terms of increasing energy,
% so |F=1,mF=1> is 1, |F=2,mF=0> is 2, etc 
%
op.initPop(1:3) = 1;
op.integrate(1e-7,50e-6);    %Integrate using time step (first argument) up to a given time (second argument) assuming constant fields

%% Plot
%
% Plots the ground state populations
%
figure(1);clf;
op.plotPopulations('ground');
hold on
P = op.getPopulations('ground');
plot(op.t,sum(P(1:3,:),1)./sum(P,1),'b--');
plot(op.t,sum(P(4:8,:),1)./sum(P,1),'r--');
h = gcf;
str = h.Children(1).String;
str{end - 1} = 'F = 1';str{end} = 'F = 2';
h.Children(1).String = str;
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
basis = 'spherical';
if strcmpi(basis,'spherical')
    str = {'-1','0','+1'};
else
    str = {'x','y','z'};
end
P = op.getPolarization(basis,'field');
P = P/op.laser1.field;
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
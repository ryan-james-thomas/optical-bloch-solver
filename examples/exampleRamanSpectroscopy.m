%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to see
%   Raman spectroscopy line shapes for square pulses. This example shows
%   Rabi flopping between the clock states of Rb-87. This will fit the
%   resulting line shape to get the Rabi frequeny and the offset due to the
%   AC Stark shift
%
op = opticalSystem('Rb87','D2');
one_photon_detuning = -1167.8679e6;
op.laser1.setGaussBeam(50e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([2,0],[3,0],one_photon_detuning);
op.laser2.setGaussBeam(50e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([1,0],[3,0],one_photon_detuning);

th = 0;
ph = 0;
op.setMagneticField(10e-3,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initPop(2) = 1;

tmp = 2*pi*1e2*ones(8,8)*0;
tmp = tmp - diag(diag(tmp));
op.decay(1:8,1:8) = tmp;


%% Solve for each detuning
f2 = 10*(-50:2:50)*1e3;
tau = 2e-6;
P = zeros(op.transition.ground.numStates,numel(f2));
Nph = zeros(numel(f2),1);
for nn = 1:numel(f2)
    op.laser2.detuning = op.laser1.detuning + f2(nn);
    op.calcBareH.integrate(tau,tau);
    tmp = op.getPopulations('ground');
    P(:,nn) = tmp(:,end);
    Nph(nn) = trapz(op.t,sum(op.getScatteringRates,1));
end

%% Plot
figure(1);clf;
ax = gca;
grid on;
set(ax,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
plot(ax,f2/1e3,P,'linewidth',0.5,'marker','.');
str = op.getPopLegend('ground');
legend(ax,str);
xlabel('Frequency [kHz]');
ylabel('Population');

%% Fit
[~,idx] = max(op.initPop);
nlf = nonlinfit(f2/1e6,P(idx,:));
nlf.useErr = false;
nlf.setFitFunc(@(A,R,x0,x) A*(1 - 4*R.^2./(4*R.^2+(x-x0).^2).*sin(sqrt(4*R.^2+(x-x0).^2)*2*pi*tau*1e6/2).^2));
[m,idx] = min(nlf.y);
Rguess = asin(sqrt(1 - m))/(2*pi*tau)*1e-6;
nlf.bounds2('A',[0.9,1,0.99],'R',[0,0.5,Rguess],'x0',[min(nlf.x),max(nlf.x),nlf.x(idx)]);
nlf.fit;
hold(ax,'on');
plot(ax,nlf.x*1e3,nlf.f(nlf.x),'k-');
hold(ax,'off');
str{9} = 'Fit';
legend(ax,str);
fprintf(1,'Rabi Freq = %.3e kHz, Center = %.3e kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);
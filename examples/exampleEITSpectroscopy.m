%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to solve
%   for steady-state spectra
%
op = opticalSystem('Rb87','D1');
op.laser1.setIntensity(1e-2*16)...
    .setPolarization([1,0,0],'linear')...
    .setStates([2,0],[2,0],1e9);
op.laser2.setIntensity(1600)...
    .setPolarization([0,1,0],'linear')...
    .setStates([1,0],[2,0],1e9);

th = 0;ph = 0;
op.setMagneticField(1,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);

D = 5e6*linspace(-1,1,500);
Pg = zeros(op.transition.ground.numStates,numel(D));
Pe = zeros(op.transition.excited.numStates,numel(D));
susc = zeros(3,numel(D));
for nn = 1:numel(D)
    op.laser1.detuning = 1e9 + 0.13e6 + D(nn);
    op.refresh.solveSteadyState;
    Pg(:,nn) = op.getPopulations('ground');
    Pe(:,nn) = op.getPopulations('excited');
    susc(:,nn) = op.getPolarization('linear');
end

%%
figure(1);clf;
ax = gca;
grid on;
set(gca,'NextPlot','ReplaceChildren','LineStyleOrder',{'-','--','.-'})
plot(ax,D/1e6,Pg,'linewidth',1.5);
legend(op.getPopLegend('ground'));
xlabel('Frequency [MHz]');
ylabel('Population');

figure(2);clf;
subplot(1,2,1);
plot(D/1e6,real(susc),'-','linewidth',2);
xlabel('Frequency [MHz]');
ylabel('Absorbance');
% legend('-1','0','+1');
legend('x','y','z');

subplot(1,2,2);
plot(D/1e6,imag(susc),'-','linewidth',2);
xlabel('Frequency [MHz]');
ylabel('Reactance');
% legend('-1','0','+1');
legend('x','y','z');


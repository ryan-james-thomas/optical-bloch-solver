%% Example script
%
%   This demonstrates how to use the optical-bloch-solver classes to show
%   some Raman Rabi flopping
%
op = opticalSystem('Rb87','D2');
op.laser1.setGaussBeam(750e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([2,0],[3,0],6.834e9 - 1.320317e9);
op.laser2.setGaussBeam(750e-3,10e-3)...
    .setPolarization([0,0,1],'spherical')...
    .setStates([1,0],[3,0],6.834e9 - 1.320317e9 + 6.882e3);

th = 0;
ph = 0;
op.setMagneticField(10e-3,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initPop(2) = 1;
%
% This bit of code attempts to calculate optical pumping rates and
% therefore decoherence associated with each laser addressing the other
% ground state manifold
%
R1 = op.getOffResonantPumping(op.laser1);
R2 = op.getOffResonantPumping(op.laser2);
R = R1 + R2;
op.makeTotalLindblad;
for gg = 1:op.transition.ground.numStates
    for ee = 1:op.transition.excited.numStates
        eShift = ee + op.transition.ground.numStates;
        op.decay(eShift,gg) = R(eShift,gg);
        op.lindblad = op.lindblad + op.makeLindblad(eShift,gg);
    end
end

tmp = 2*pi*1e2*ones(8,8)*0;
tmp = tmp - diag(diag(tmp));
op.decay(1:8,1:8) = tmp;

op.integrate(0.1e-6,10e-6);

t0 = op.t;
%%
op.t = t0*1e6;
figure(4);clf;
subplot(2,1,1);
op.plotPopulations('ground');
xlabel('Time [us]');
xlim([0,Inf]);
subplot(2,1,2);
P = op.getPopulations('ground');
plot(op.t,sum(P(1:3,:),1),'k-','linewidth',2);
hold on;
plot(op.t,sum(P(4:8,:),1),'r-','linewidth',2);
xlabel('Time [us]');
ylabel('Populations');
legend('F = 1','F = 2');
xlim([0,Inf]);

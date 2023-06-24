%% Example script
%
%   This calculates the effective OD for a sample of atoms
%
dt = 1e-6;
T_repump = 200e-6;
I_repump = 16;
T_imaging = 40e-6;
P_imaging = 10e-6;
w_imaging = 2.3e-3;
op = opticalSystem('Rb87','D2');
th = 0*pi/180;ph = 0;
B = 0.1;

% D = 1e6*linspace(-50,50,51);
D = 1e6*(-12:12);

%% Pump from |1,-1> to F = 2 manifold
op.laser1.setIntensity(I_repump)...
    .setPolarization([1,1,1],'linear')...
    .setStates([1,-1],[2,-2],0);

op.setMagneticField(B,[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]);
op.initPop(3) = 1;
op.integrate(T_repump,T_repump);

P = op.getPopulations('ground');
Pf = P(:,end);

%% Calculate for each detuning
op.laser1.setStates([2,0],[3,0],0)...
    .setPolarization([0,0,1],'spherical')...
    .setGaussBeam(P_imaging,w_imaging);

op.initPop(1:8) = Pf;
% op.initPop(1:8) = 0;op.initPop(4) = 1;
polarization = zeros(numel(D),3);
for nn = 1:numel(D)
    op.laser1.detuning = D(nn);
    op.refresh.integrate(dt,T_imaging);
    pol = op.getPolarization('spherical');
    polarization(nn,:) = trapz(op.t,pol,2);
    
%     figure(1);clf;
%     P = op.getPopulations('ground');
%     plot(op.t*1e6,P(4:8,:),'-','linewidth',2);
%     yyaxis right;
%     plot(op.t*1e6,-real(pol(3,:))/op.laser1.field,'-','linewidth',2);
%     pause(1e-3);
end

%% Plot
cloud_size = sqrt(const.kb*1.5e-6/const.mRb)*19e-3;
% cloud_size = 5e-3;
Natoms = 1e6;
n0 = Natoms/(2*pi*cloud_size^2);

OD = -real(polarization(:,3))/op.laser1.field*n0/T_imaging;

nlf = nonlinfit(D/1e6,OD,1e-2);
nlf.setFitFunc(@(A,x0,w,y0,x) y0 + A./(1 + 4*(x - x0).^2/w^2));
nlf.bounds2('A',[0,2,max(nlf.y)],'x0',[min(nlf.x),max(nlf.x),0],'w',[0,20,6],'y0',[0,0.05,0.001]);
nlf.fit

figure(3);clf;
% plot(D/1e6,OD,'o-');
nlf.plot('plotresiduals',false);
plot_format('Detuning [MHz]','Optical depth','',10);

disp('Creating simulated data with homogeneous background conductivity');

vacuum_perm = 8.85e-12;
perm = 4*vacuum_perm;
cond = 2e-3;
d1 = 5e-3;
d2 = 2e-3;
omega = 3.2e5;
coupl1 = perm*omega/d1;
coupl2 = perm*omega/d2;
couplw = coupl1*10;

err = [1e-3 0 0 0];

coupl0 = 1e-9;
cond1 = cond;
cond2 = cond;

load('meshfiles/DLSS_irto_densest.mat');%Load mesh variables
load('meshfiles/DLSS_irto_inv.mat');
%load('meshfiles/DLSS_irto_dense.mat');%Load mesh variables
%load('meshfiles/DLSS_irto_inv.mat');
ginvsimu = ginv;
Hinvsimu = Hinv;
coupling = GenerateEllipseInGrid(ginvsimu, Hinvsimu, coupl1, coupl2, 5e-2, 5e-2, 10e-2, 15e-2, 2e-2, 1);%Generate a single blob as a target
z = 1e-6*ones(length(elfaces),1);
offset = [1 0];
ng1 = length(g1);
ng0 = length(g0);
H1 = H0;
H2 = H0+ng1;
for ie = 1:length(elfaces)
    if els(ie) == 1
        H1 = [H1; Hel{ie}];
        elfaces{ie} = elfaces{ie}+ng1;
    elseif els(ie) == 2
        H2 = [H2; Hel{ie}+ng1];
    end
end
eln = zeros(size(H1,1) + size(H2,1),1);
counter1 = size(H0,1)+1;
counter2 = size(H1,1) + size(H0,1)+1;
for ie = 1:length(Hel)
    tempn = size(Hel{ie},1);
    if els(ie) == 1
        eln(counter1:counter1+tempn-1) = ie;
        counter1 = counter1 + tempn;
    elseif els(ie) == 2
        eln(counter2:counter2+tempn-1) = ie;
        counter2 = counter2 + tempn;
    end
end
ci = zeros(ng1 + length(g2),1);
ci(1:ng0) = (1:ng0)' + ng1;
ci(ng1+1:ng1+ng0) = (1:ng0)';
fg = [g1; g2+offset];
fH = [H1; H2];
fmsimu = ForwardMeshDLSS(fg, fH, elfaces, ng1, ci, eln, offset);
imesh.g = ginvsimu;
imesh.H = Hinvsimu;
load('meshfiles/DLSS_irto_Layer1.mat');
sinvg = g;
sinvH = H;
nsinvg = size(g,1);
load('meshfiles/DLSS_irto_Layer2.mat');
simesh.g = [sinvg; g+offset];
simesh.H = [sinvH; H+nsinvg];
fmsimu.SetInverseMesh(simesh);
sigmasimu = [cond1*ones(nsinvg,1); cond2*ones(length(simesh.g)-nsinvg,1)];

pl = Plotter(simesh.g, simesh.H);
pl.plot(sigmasimu);
axis equal;
axis off;
ccc = colorbar;
ccc.FontSize = 16;
saveas(pl.fig, 'simulations/target_W1_cond.png');

plc = Plotter(ginvsimu, Hinvsimu);
plc.plot(coupling);
axis off;
ccc = colorbar;
ccc.FontSize = 16;
saveas(plc.fig, 'simulations/target_W1_coupl.png');

[Umeas, Imeas, Umeas_i, Imeas_i, Umeas_0, Imeas_0] = Simulation_DLSS_three_phase(fmsimu, sigmasimu, coupling, z, vincl, imesh, err, mode, 0, coupl1, nsinvg);
save([simuname '_W1.mat'], 'Umeas', 'Imeas', 'sigmasimu', 'z', 'ginvsimu', 'Hinvsimu', 'Imeas_i', 'Umeas_i', 'Imeas_0', 'Umeas_0', 'coupling');
disp(['simulation run succesfully, results saved in ' simuname '_W1.mat']);
reco.sigma = sigmasimu;
reco.coupl = coupling;
reco.zeta = z;
save('simu_W1.mat', 'reco');

%==================================================================

disp('Creating simulated data with inhomogeneous background conductivity');

coupling = GenerateEllipseInGrid(ginvsimu, Hinvsimu, coupl1, coupl2, 5e-2, 5e-2, 20e-2, 10e-2, 2e-2, 1);%Generate a single blob as a target

SmoothnessSigma = PriorSmoothness(simesh.g, 0.1, (2e-1*cond)^2, cond*ones(size(simesh.g,1),1));
sigmasimu = SmoothnessSigma.DrawSamples(1);
sigmasimu(sigmasimu<1e-2*cond) = 1e-2*cond;

pl.plot(sigmasimu);
axis equal;
axis off
ccc = colorbar;
ccc.FontSize = 16;
saveas(pl.fig, 'simulations/target_W2_cond.png');

plc.plot(coupling);
axis off
ccc = colorbar;
ccc.FontSize = 16;
saveas(plc.fig, 'simulations/target_W2_coupl.png');

z = randn(2*length(elfaces),1);
z(z<1e-6) = 1e-6;
figure(); plot(z); saveas(gcf, 'simulations/target_W2_zeta.png');

[Umeas, Imeas, Umeas_i, Imeas_i, Umeas_0, Imeas_0] = Simulation_DLSS_three_phase(fmsimu, sigmasimu, coupling, z, vincl, imesh, err, mode, 0, coupl1, nsinvg);
save([simuname '_W2.mat'], 'Umeas', 'Imeas', 'sigmasimu', 'z', 'ginvsimu', 'Hinvsimu', 'Imeas_i', 'Umeas_i', 'Imeas_0', 'Umeas_0', 'coupling');
disp(['simulation run succesfully, results saved in ' simuname '_W2.mat']);

reco.sigma = sigmasimu;
reco.coupl = coupling;
reco.zeta = z;
save('simu_W2.mat', 'reco');

%=================================================================

disp('Creating simulated data of water injection');

coupling = GenerateEllipseInGrid(ginvsimu, Hinvsimu, coupl1, couplw, 10e-2, 10e-2, 15e-2, 5e-2, 5e-2, 1);%Generate a single blob as a target


SmoothnessSigma = PriorSmoothness(simesh.g, 0.1, (2e-1*cond)^2, cond*ones(size(simesh.g,1),1));
sigmasimu = SmoothnessSigma.DrawSamples(1);
sigmasimu(sigmasimu<1e-2*cond) = 1e-2*cond;
deltas1 = GenerateEllipseInGrid(simesh.g, simesh.H, 0, cond, 10e-2, 10e-2, 15e-2, 5e-2, 5e-2, 1);%Generate a single blob as a target
deltas2 = GenerateEllipseInGrid(simesh.g, simesh.H, 0, cond, 10e-2, 10e-2, 15e-2+offset(1), 5e-2+offset(2), 5e-2, 1);%Generate a single blob as a target
sigmasimu2 = sigmasimu + deltas1 + deltas2;

pl.plot(sigmasimu);
axis equal;
axis off
ccc = colorbar;
ccc.FontSize = 16;
saveas(pl.fig, 'simulations/target_water_cond.png');

pl.plot(sigmasimu2);
axis equal;
axis off
ccc = colorbar;
ccc.FontSize = 16;
saveas(pl.fig, 'simulations/target_water_cond2.png');


plc.plot(coupling);
axis off
ccc = colorbar;
ccc.FontSize = 16;
saveas(plc.fig, 'simulations/target_water_coupl.png');

z = randn(2*length(elfaces),1);
z(z<1e-6) = 1e-6;
figure(); plot(z); saveas(gcf, 'simulations/target_water_zeta.png');

[Umeas, Imeas, Umeas_i, Imeas_i, Umeas_0, Imeas_0] = Simulation_DLSS_three_phase(fmsimu, sigmasimu, coupling, z, vincl, imesh, err, mode, 0, coupl1, nsinvg, sigmasimu2);
save([simuname '_water.mat'], 'Umeas', 'Imeas', 'sigmasimu', 'z', 'ginvsimu', 'Hinvsimu', 'Imeas_i', 'Umeas_i', 'Imeas_0', 'Umeas_0', 'coupling', 'sigmasimu2');
disp(['simulation run succesfully, results saved in ' simuname '_water.mat']);

reco.sigma = sigmasimu;
reco.coupl = coupling;
reco.zeta = z;
save('simu_water.mat', 'reco');

reco.sigma = sigmasimu2;
reco.coupl = coupling;
reco.zeta = z;
save('simu_water2.mat', 'reco');



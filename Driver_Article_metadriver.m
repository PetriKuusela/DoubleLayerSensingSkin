clear all;
close all;

addpath('./OOEIT1.0');

NEW_DATA = 1;
COMPUTE_RECONSTRUCTIONS_SIMU = 1;
COMPUTE_RECONSTRUCTIONS_REAL = 1;
FORMAT_FIGURES = 1;

load_cond_data = 0;
load_epshest_data = 0;
load_sigma_data_real = 0;
load_epshest_data_real = 0;
coupl0 = 0;

simuname = 'Article_simus';

mode = 'potential';
vincl = true(32,32);
vincl = vincl(:);
cond = 2e-3;

if NEW_DATA
    Generate_data;
end

%% Load meshes for inverse problem
%Load meshes for solving the inverse problem, i.e. a 3D forward mesh and a
%2D inverse mesh.
load('meshfiles/DLSS_irto_inv.mat');
load('meshfiles/DLSS_irto.mat');
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
        %elfaces{ie} = elfaces{ie}+ng1;
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
fm = ForwardMeshDLSS([g1; g2+offset], [H1; H2], elfaces, ng1, ci, eln, offset);
felfaces = elfaces;
load('meshfiles/DLSS_irto_Layer1.mat');
sinvg = g;
sinvH = H;
nsinvg = size(g,1);
load('meshfiles/DLSS_irto_Layer2.mat');
simesh.g = [sinvg; g+offset];
simesh.H = [sinvH; H+nsinvg];
disp('meshfiles loaded');



if COMPUTE_RECONSTRUCTIONS_SIMU

    efc = 1e-4;
    load([simuname '_W1.mat'])
    Driver_estimate_sigma_auto;
    save('simulations/W1.mat', 'reco', 'InvSolver');
    saveas(InvSolver.Plotter.fig, 'simulations/W1_cond.png');
    saveas(InvSolver.Plotter.fig2, 'simulations/W1_coupl.png');
    clear InvSolver reco;
    close all;
    disp('Completed W1');

    load([simuname '_W2.mat'])
    Driver_estimate_sigma_auto;
    save('simulations/W2.mat', 'reco', 'InvSolver');
    saveas(InvSolver.Plotter.fig, 'simulations/W2_cond.png');
    saveas(InvSolver.Plotter.fig2, 'simulations/W2_coupl.png');
    clear InvSolver reco;
    close all;

    disp('Completed W2');

    efc = 1e-2;
    load([simuname '_water.mat'])
    Driver_estimate_sigma_auto;
    save('simulations/water.mat', 'reco', 'InvSolver');
    saveas(InvSolver.Plotter.fig, 'simulations/water_cond.png');
    saveas(InvSolver.Plotter.fig2, 'simulations/water_coupl.png');
    clear InvSolver reco;
    close all;

    disp('Completed Water');
end

if COMPUTE_RECONSTRUCTIONS_REAL
    datanames = {'Irto_DLSS2/Layers_ref', 'Irto_DLSS2/Layers_M1_ref2', 'Irto_DLSS2/Layers_M1_W2'; 
                 'Irto_DLSS2/Layers_ref', 'Irto_DLSS2/Layers_M1_ref3', 'Irto_DLSS2/Layers_M1_W3'; 
                 'Irto_DLSS2/Layers_ref', 'Irto_DLSS2/Layers_M4_ref3', 'Irto_DLSS2/Layers_M4_water4'};

    for iim = 1:size(datanames,1)
        Driver_realdata_auto;
        save(['realdata/r' num2str(iim) '.mat'], 'reco', 'InvSolver');
        saveas(InvSolver.Plotter.fig, ['realdata/r' num2str(iim) '_cond.png']);
        saveas(InvSolver.Plotter.fig2, ['realdata/r' num2str(iim) '_coupl.png']);
        clear InvSolver reco;
        close all;
        disp(['Completed realdata n. ' num2str(iim)]);
    end
end

if FORMAT_FIGURES

    %Very unnecessary, but still, reload the mesh variables with the 2
    %layers closer together for plotting purposes
    load('meshfiles/DLSS_irto_inv.mat');
    load('meshfiles/DLSS_irto.mat');
    offset = [0.3 0];
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
            %elfaces{ie} = elfaces{ie}+ng1;
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
    fm = ForwardMeshDLSS([g1; g2+offset], [H1; H2], elfaces, ng1, ci, eln, offset);
    felfaces = elfaces;
    load('meshfiles/DLSS_irto_Layer1.mat');
    sinvg = g;
    sinvH = H;
    nsinvg = size(g,1);
    load('meshfiles/DLSS_irto_Layer2.mat');
    simesh.g = [sinvg; g+offset];
    simesh.H = [sinvH; H+nsinvg];

    pl = Plotter_estsig(simesh.g, simesh.H, ginv, Hinv);

    files = {'simu_W1.mat'; 'simulations/W1.mat'; 'simu_W2.mat'; 'simulations/W2.mat'};
    cbname = 'simulations/simu_cb.png';
    scbname = ['simulations/simu_cb_s.png'];
    Figureformatter_auto;

    files = {'simu_water.mat'; 'simu_water2.mat'; 'simulations/water.mat'};
    cbname = 'simulations/simu_cb_w.png';
    scbname = ['simulations/simu_cb_w_s.png'];
    Figureformatter_auto;

    files = {'realdata/r1.mat', 'realdata/r2.mat'};
    cbname = 'realdata/cb.png';
    scbname = 'realdata/cb_s.png';
    Figureformatter_auto;

    files = {'realdata/r3.mat'};
    cbname = 'realdata/cb_w.png';
    scbname = 'realdata/cb_w_s.png';
    Figureformatter_auto;


    % files = {'simulations/W1.mat'; 'simulations/W2.mat'};
    % cbname = 'simulations/cb.png';
    % scbname = 'simulations/cb_s.png';
    % Figureformatter_auto;
    % 
    % files = {'simulations/water.mat'};
    % cbname = 'simulations/cb_w.png';
    % scbname = 'simulations/cb_w_s.png';
    % Figureformatter_auto;
    % 
    % files = {'realdata_better/r2.mat', 'realdata_better/r3.mat'};
    % cbname = 'realdata_better/cb.png';
    % scbname = 'realdata_better/cb_s.png';
    % Figureformatter_auto;
    % 
    % files = {'realdata_better/r11.mat'};
    % cbname = 'realdata_better/cb_w.png';
    % scbname = 'realdata_better/cb_w_s.png';
    % Figureformatter_auto;
end


%% Loading data
disp('Loading data');

l1name = datanames{iim, 1};
refname = datanames{iim, 2};
dataname = datanames{iim, 3};

dataname = [dataname '/ect_1/data1'];
refname = [refname '/ect_1/data1'];
l1name = [l1name '/ect_1/data1'];

coupl = 0.35;
cond = 2e-3;
mode = 'potential';

vincl = true(32,32);
vincl = vincl(:);
nel16 = 32;

vincl16 = true(nel16,nel16);
vincl16 = vincl16(:);

disp(['Loading data from file ' dataname ', and reference data from file ' refname '.']);
addpath('../OOEIT1.0/Rocsole');
[Umeas_si1, Imeas_si1] = Dataloader(l1name, vincl16, nel16, nel16, 10);
[Umeas_i, Imeas_i] = Dataloader(refname, vincl, 32, 32, 10);
[Umeas, Imeas] = Dataloader(dataname, vincl, 32, 32, 10);

Umeas_si = zeros(32,32);
Imeas_si = zeros(32,32);
Imeas_si1 = reshape(Imeas_si1, nel16, nel16);
Umeas_si(2:2:32,2:2:32) = Umeas_si1(1:16,1:16);
Umeas_si(1:2:31,1:2:31) = Umeas_si1(17:32,17:32);
Imeas_si(2:2:32,2:2:32) = Imeas_si1(1:16,1:16);
Imeas_si(1:2:31,1:2:31) = Imeas_si1(17:32,17:32);
Umeas_si = Umeas_si(:);
Imeas_si = Imeas_si(:);

Imeas_i = reshape(Imeas_i, 32, 32);
Umeas_ti(2:2:32,2:2:32) = Umeas_i(1:16,1:16);
Umeas_ti(1:2:31,1:2:31) = Umeas_i(17:32,17:32);
Umeas_ti(2:2:32,1:2:31) = Umeas_i(1:16,17:32);
Umeas_ti(1:2:31,2:2:32) = Umeas_i(17:32,1:16);
Imeas_ti(2:2:32,2:2:32) = Imeas_i(1:16,1:16);
Imeas_ti(1:2:31,1:2:31) = Imeas_i(17:32,17:32);
Imeas_ti(2:2:32,1:2:31) = Imeas_i(1:16,17:32);
Imeas_ti(1:2:31,2:2:32) = Imeas_i(17:32,1:16);
Umeas_i = Umeas_ti(:);
Imeas_i = Imeas_ti(:);

Imeas = reshape(Imeas, 32, 32);
Umeas_t(2:2:32,2:2:32) = Umeas(1:16,1:16);
Umeas_t(1:2:31,1:2:31) = Umeas(17:32,17:32);
Umeas_t(2:2:32,1:2:31) = Umeas(1:16,17:32);
Umeas_t(1:2:31,2:2:32) = Umeas(17:32,1:16);
Imeas_t(2:2:32,2:2:32) = Imeas(1:16,1:16);
Imeas_t(1:2:31,1:2:31) = Imeas(17:32,17:32);
Imeas_t(2:2:32,1:2:31) = Imeas(1:16,17:32);
Imeas_t(1:2:31,2:2:32) = Imeas(17:32,1:16);
Umeas = Umeas_t(:);
Imeas = Imeas_t(:);


%% Setting up the inverse solver

%Set up the forward problem solver:
imesh.g = ginv;
imesh.H = Hinv;
fm.SetInverseMeshCoupling(imesh);
fm.SetInverseMesh(simesh);
fname = 'autosave';
if load_sigma_data_real
    load([fname '_cond.mat']);
    disp(['loaded epsilon correction from file ' fname '_cond.mat']);
else
    solver = EITFEM_DLSS_multiest(fm);
    solver.sigmamin = 1e-5;
    solver.zeta = 1e-6*ones(length(felfaces),1);
    solver.mode = mode;
    solver.sigma = cond*ones(size(simesh.g,1),1);
    solver.coupl = coupl0*ones(size(ginv,1),1);
    solver.vincl = vincl;
    if strcmp(mode, 'current')
        solver.IDmeas = Imeas_si;
        solver.Dmeas = Umeas_si;
    elseif strcmp(mode, 'potential')
        solver.IDmeas = Umeas_si;
        solver.Dmeas = Imeas_si;
    end
    solver.SetInvGamma(1e-4, 3e-2);%(1e-4, 3e-2);
    disp('solving epsilon correction and conductivity estimate');

    %Set up the smoothness prior:
    gstruct.coupl = ginv;
    gstruct.sigma = simesh.g;
    meanstruct.coupl = zeros(size(ginv,1),1);
    meanstruct.sigma = cond*ones(size(simesh.g,1),1);
    varstruct.coupl = 1e-3;
    varstruct.sigma = (1e-2*cond)^2;%1e-1;
    SmoothPrior = PriorSmoothness_estsig(gstruct, 0.05, varstruct, meanstruct);
    disp('Smoothness prior for conductivity estimation set up');

    %Set up the positivity prior:
    mastruct.sigma = 1e-5;
    mastruct.coupl = 1e-6;
    mastruct.zeta = 1e-5;
    va0struct.sigma = 1e-1;
    va0struct.coupl = 1e-3;
    va0struct.zeta = 1e-2;
    PosiPrior = PriorPositivityParabolic_estsig(mastruct, va0struct);
    disp('Positivity prior for conductivity estimation set up');

    plotter = Plotter_estsig(simesh.g, simesh.H, ginv, Hinv);

    resobj = cell(3, 1);
    resobj{1} = solver;
    resobj{2} = SmoothPrior;
    resobj{3} = PosiPrior;
    InvSolver = SolverLinesearch(resobj);
    InvSolver.maxIter = 10;
    InvSolver.Plotter = plotter;
    InvSolver.showSplitVals = 1;
    
    %Make the initial guess and start solving!
    %sigmainitial.coupl = coupl0*ones(size(ginv,1),1);
    sigmainitial.sigma = [cond*ones(nsinvg,1); cond*ones(size(simesh.g,1)-nsinvg, 1)];
    zetainitial = 1e-2*ones(32,1);
    sigesti = DLSS_est(sigmainitial.sigma, [], zetainitial);
    disp('Beginning to estimate the conductivity with 0 coupling.');
    condreco = InvSolver.solve(sigesti);
    disp('Conductivity reconstruction computed!');

    %vals = solver.SolveForwardVec(condreco);
    %eps = vals - Imeas_i;
    save([fname '_cond.mat'], 'condreco');
end

if load_epshest_data_real
    load([fname '_epshest.mat']);
    disp(['loaded epsilon correction from file ' fname '_epshest.mat']);
else

    solver = EITFEM_DLSS_multiest(fm);
    solver.sigmamin = 1e-5;
    solver.zeta = condreco.zeta;%1e-6*ones(length(felfaces),1);
    solver.sigma = condreco.sigma;
    solver.vincl = vincl;
    if strcmp(mode, 'current')
        solver.IDmeas = Imeas_i;
        solver.Dmeas = Umeas_i;
    elseif strcmp(mode, 'potential')
        solver.IDmeas = Umeas_i;
        solver.Dmeas = Imeas_i;
    end
    solver.SetInvGamma(1e-4, 3e-2);%(1e-4, 3e-2);
    sigesti = DLSS_est([], ones(size(ginv,1),1), []);
    [eps, couplhest] = SolveEpsilonCorrection(solver, sigesti);

    %vals = solver.SolveForwardVec(condreco.sigma, condreco.coupl);
    %eps = vals - Imeas_si;
    save([fname '_epshest.mat'], 'eps','couplhest');
end

cond = mean(condreco.sigma);

solver = EITFEM_DLSS_multiest(fm);
solver.sigmamin = 1e-5;
solver.zeta = condreco.zeta;%1e-6*ones(length(felfaces),1);
solver.mode = mode;
solver.sigma = condreco.sigma;
solver.coupl = couplhest*ones(size(imesh.g,1),1);
imesh.g = ginv;
imesh.H = Hinv;
fm.SetInverseMeshCoupling(imesh);
fm.SetInverseMesh(simesh);
solver.vincl = vincl;
solver.eps = eps;

if strcmp(mode, 'current')
    solver.IDmeas = Imeas;
    solver.Dmeas = Umeas;
elseif strcmp(mode, 'potential')
    solver.IDmeas = Umeas;
    solver.Dmeas = Imeas;
end
solver.SetInvGamma(1e-4, 3e-2);
disp('Forward problem solver set up')

%Set up the Total Variation prior:
ec.coupl = couplhest;%expected change
ec.sigma = cond;
ef.coupl = 1;
ef.sigma = 1e-5;
TVPrior = PriorTotalVariation_estsig(simesh.g, simesh.H, ginv, Hinv, ec, ef, 1e-5);
disp('TV prior set up');

%Set up the smoothness prior:
gstruct.coupl = ginv;
gstruct.sigma = simesh.g;
meanstruct.coupl = couplhest*ones(size(ginv,1),1);
meanstruct.sigma = condreco.sigma;
varstruct.coupl = 1e6;%couplhest^2;%1e6;
varstruct.sigma = (1e-2*cond)^2;%(1e-1*cond)^2;
SmoothPrior = PriorSmoothness_estsig(gstruct, 0.05, varstruct, meanstruct);
disp('Smoothness prior set up');

%Set up the positivity prior:
mastruct.sigma = 1e-4;
mastruct.coupl = 1e-2;
mastruct.zeta = 1e-6;
va0struct.sigma = 1e1;
va0struct.coupl = 1e2;
va0struct.zeta = 1e-1;
PosiPrior = PriorPositivityParabolic_estsig(mastruct, va0struct);
disp('Positivity prior set up');

ev = condreco.zeta;
strength = 1e-2;
ZetaPrior = PriorZeta_estsig(ev, strength);


%Set up the plotter:
plotter = Plotter_estsig(simesh.g, simesh.H, ginv, Hinv);

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(4, 1);
resobj{1} = solver;
resobj{2} = SmoothPrior;
resobj{3} = TVPrior;
resobj{4} = PosiPrior;
%resobj{5} = ZetaPrior;
InvSolver = SolverLinesearch(resobj);
InvSolver.maxIter = 150;
InvSolver.Plotter = plotter;
InvSolver.showSplitVals = 1;
InvSolver.erel = 1e-6;

%Make the initial guess and start solving!
sigmainitial.coupl = couplhest*ones(size(ginv,1),1);
sigmainitial.sigma = condreco.sigma;
sigmainitial.zeta = condreco.zeta;
sigesti = DLSS_est(sigmainitial.sigma, sigmainitial.coupl, []);
disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(sigesti);        






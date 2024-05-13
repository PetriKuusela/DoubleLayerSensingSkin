%% Setting up the inverse solver

%Set up the forward problem solver:
imesh.g = ginv;
imesh.H = Hinv;
fm.SetInverseMeshCoupling(imesh);
fm.SetInverseMesh(simesh);
fname = 'autosave';
if load_cond_data
    load([fname '_cond.mat']);
    disp(['loaded conductivity estimates from file ' fname '_cond.mat']);
else
    solver = EITFEM_DLSS_multiest(fm);
    solver.sigmamin = 1e-9;
    solver.zeta = 1e-6*ones(length(felfaces),1);
    solver.mode = mode;
    solver.sigma = cond*ones(size(simesh.g,1),1);
    solver.coupl = zeros(size(ginv,1),1);
    solver.vincl = vincl;
    if strcmp(mode, 'current')
        solver.IDmeas = Imeas_0;
        solver.Dmeas = Umeas_0;
    elseif strcmp(mode, 'potential')
        solver.IDmeas = Umeas_0;
        solver.Dmeas = Imeas_0;
    end
    solver.SetInvGamma(1e-4, 3e-2);%(1e-4, 3e-2);

    %Set up the smoothness prior:
    gstruct.coupl = ginv;
    gstruct.sigma = simesh.g;
    meanstruct.coupl = zeros(size(ginv,1),1);
    meanstruct.sigma = cond*ones(size(simesh.g,1),1);
    varstruct.coupl = 1e-3;
    varstruct.sigma = (5e-1*cond)^2;%1e-1;
    SmoothPrior = PriorSmoothness_estsig(gstruct, 0.2, varstruct, meanstruct);
    disp('Smoothness prior for conductivity estimation set up');

    %Set up the positivity prior:
    mastruct.sigma = 1e-5;
    mastruct.coupl = 1e-6;
    mastruct.zeta = 1e-7;
    va0struct.sigma = 1e-3;
    va0struct.coupl = 1e-3;
    va0struct.zeta = 1e-1;
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
    %sigmainitial.coupl = zeros(size(ginv,1),1);
    sigmainitial.sigma = [cond*ones(nsinvg,1); cond*ones(size(simesh.g,1)-nsinvg, 1)];
    zetainitial = 1e-6*ones(32,1);
    sigesti = DLSS_est(sigmainitial.sigma, [], zetainitial);
    disp('Beginning to estimate the conductivity with 0 coupling.');
    condreco = InvSolver.solve(sigesti);
    disp('Conductivity reconstruction computed!');

    %vals = solver.SolveForwardVec(condreco);
    %eps = vals - Imeas_i;
    save([fname '_cond.mat'], 'condreco');
end


if load_epshest_data
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

    save([fname '_epshest.mat'], 'eps','couplhest');
end

condhest = mean(condreco.sigma);

solver = EITFEM_DLSS_multiest(fm);
solver.sigmamin = 1e-9;
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
solver.SetInvGamma(1e-4, 3e-2);%(1e-4, 3e-2);
disp('Forward problem solver set up')

%Set up the Total Variation prior:
ec.coupl = couplhest;%expected change
ec.sigma = condhest;
ef.coupl = efc;%1e-4 for weights, 1e-2 for water;
ef.sigma = 1e-7;
TVPrior = PriorTotalVariation_estsig(simesh.g, simesh.H, ginv, Hinv, ec, ef);
disp('TV prior set up');

%Set up the smoothness prior:
gstruct.coupl = ginv;
gstruct.sigma = simesh.g;
meanstruct.coupl = couplhest*ones(size(ginv,1),1);
meanstruct.sigma = condreco.sigma;
varstruct.coupl = 1e9;
varstruct.sigma = (3e-1*cond)^2;%1e-1;
SmoothPrior = PriorSmoothness_estsig(gstruct, 0.2, varstruct, meanstruct);
disp('Smoothness prior set up');

%Set up the positivity prior:
mastruct.sigma = 1e-5;
mastruct.coupl = 1e-6;
mastruct.zeta = 1e-10;
va0struct.sigma = 1e-3;
va0struct.coupl = 1e-3;
va0struct.zeta = 1e-1;
PosiPrior = PriorPositivityParabolic_estsig(mastruct, va0struct);
disp('Positivity prior set up');

%Set up the plotter:
plotter = Plotter_estsig(simesh.g, simesh.H, ginv, Hinv);

%Finally, set up the inverse problem solver, in this case GN with
%linesearch:
resobj = cell(4, 1);
resobj{1} = solver;
resobj{2} = SmoothPrior;
resobj{3} = TVPrior;
resobj{4} = PosiPrior;
InvSolver = SolverLinesearch(resobj);
InvSolver.maxIter = 150;
InvSolver.Plotter = plotter;
InvSolver.showSplitVals = 1;
InvSolver.erel = 1e-6;
InvSolver.estep = 1e-12;

%Make the initial guess and start solving!
sigmainitial.coupl = couplhest*ones(size(ginv,1),1);
sigmainitial.sigma = condreco.sigma;
sigmainitial.zeta = condreco.zeta;%1e-6*ones(32,1);
sigesti = DLSS_est(sigmainitial.sigma, sigmainitial.coupl, []);
disp('All set! Beginning to solve the inverse problem.')
reco = InvSolver.solve(sigesti);            






function [eps, hest] = SolveEpsilonCorrection_DLSS_onlycoupl(solver, sigmang1, hinit)
    %This function fits a homogeneous estimate (conductivity and 
    %coupling) to the data of EITFEM_DLSS-solver, and
    %then proceeds to calculate the epsilon correction (i.e. the
    %difference between the results of FEM and measurements).
    %
    %This function utilizes SolverLinesearch-class to solve the 2D
    %optimization problem.
    %
    %Author: Petri Kuusela 16.11.2021
    
    if nargin < 3
        hinit = 1;
    end
    
    
    hw = HomogeneousWrapper(solver);
    PosiPrior = PriorPositivityParabolic(1e-6, 1e3);
    resobj = cell(2,1);
    resobj{1} = hw;
    resobj{2} = PosiPrior;
    InvSolver = SolverLinesearch(resobj);
    InvSolver.maxIter = 4;
    InvSolver.plotIterations = 0;
    InvSolver.showSplitVals = 1;
    
    if 0 == 1%strcmp(solver.mode, 'potential')
        solver.mode = 'current';
        solver.Imeas = -solver.Imeas;
        hinit = InvSolver.solve(hinit);%if the initial guess is really far, the solver struggles with mode 'potential'
        InvSolver.maxIter = 3;
        solver.mode = 'potential';
        solver.Imeas = -solver.Imeas;
    end
    
    
    hest = InvSolver.solve(hinit);
    
    elval = solver.SolveForwardVec(solver.sigma, hest*ones(hw.ngc,1));
    eps = solver.Dmeas-elval;
    
end
classdef HomogeneousWrapper < handle
%This is a wrapper class for EITFEM_DLSS class used for computing the
%homogeneous estimates. It changes the input from being full coupling to
%just one value, the homogeneous coupling
    
    
    properties
        solver  %This is the EITFEM_DLSS object to be wrapped
        ng      %The number of nodes in the inverse mesh
        ngc
    end
    
    methods
        
        function obj = HomogeneousWrapper(solver)
            %Class constructor. Input argument is the solver to be wrapped.
            obj.solver = solver;
            if numel(solver.fmesh.itof) == 1%This means that the inverse and forward mesh are the same
                obj.ng = solver.fmesh.ng;
            else%Otherwise we can get the inverse size mesh from solver.P1st
                obj.ng = size(solver.fmesh.itof,2);
            end
            if numel(solver.fmesh.itofc) == 1%This means that the inverse and forward mesh are the same
                warning('Check this before using!');
                obj.ngc = solver.fmesh.ng;
            else%Otherwise we can get the inverse size mesh from solver.P1st
                obj.ngc = size(solver.fmesh.itofc,2);
            end
        end
        function res = OptimizationFunction(self,hest)
            %self.solver.sigma = hest(2)*ones(self.ng,1);
            res = self.solver.OptimizationFunction(hest*ones(self.ngc,1));
        end
        function [Hess, grad] = GetHessAndGrad(self,hest)
            %self.solver.sigma = hest(2)*ones(self.ng,1);
            %self.solver.couplforJ = hest(1)*ones(self.ngc,1);
            [Hess, grad] = self.solver.GetHessAndGrad(hest(1)*ones(self.ngc,1));
            %[Hess2, grad2, Hessmix] = self.solver.GetHessAndGradSigma(hest(1)*ones(self.ngc,1), hest(2)*ones(self.ng,1));
            %grad = [sum(grad); sum(grad2)];
            grad = sum(grad);
            Hess = sum(sum(Hess));
        end
        function Plot(self,hest)
            %This can be used to plot the fitting of the data during the
            %solving of the homogeneous estimate.
            %self.solver.sigma = hest(2)*ones(self.ng,1);
            self.solver.Plot(hest*ones(self.ngc,1));
        end
        
        
    end
    
end
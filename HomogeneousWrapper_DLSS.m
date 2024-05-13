classdef HomogeneousWrapper_DLSS < handle
%This is a wrapper class for EITFEM_DLSS class used for computing the
%homogeneous estimates. It changes the input from being full coupling to
%two values: constant coupling and conductivity.
    
    
    properties
        solver  %This is the EITFEM_DLSS object to be wrapped
        ng      %The number of nodes in the inverse mesh
        ngc
        sng1     %the number of nodes in first layer of sigma mesh
    end
    
    methods
        
        function obj = HomogeneousWrapper_DLSS(solver, sng1)
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
            obj.sng1 = sng1;
        end
        function res = OptimizationFunction(self,hest)
            self.solver.sigma = [hest(2)*ones(self.sng1,1); hest(3)*ones(self.ng-self.sng1,1)];
            res = self.solver.OptimizationFunction(hest(1)*ones(self.ngc,1));
        end
        function [Hess, grad] = GetHessAndGrad(self,hest)
            self.solver.sigma = [hest(2)*ones(self.sng1,1); hest(3)*ones(self.ng-self.sng1,1)];
            self.solver.couplforJ = hest(1)*ones(self.ngc,1);
            [Hess, grad] = self.solver.GetHessAndGrad(hest(1)*ones(self.ngc,1));
            [Hess2, grad2, Hessmix] = self.solver.GetHessAndGradSigma(hest(1)*ones(self.ngc,1), self.solver.sigma);
            grad = [sum(grad); sum(grad2(1:self.sng1)); sum(grad2(self.sng1+1:end))];
            Hess = [sum(sum(Hess)) sum(sum(Hessmix(1:self.sng1,:))) sum(sum(Hessmix(self.sng1+1:end,:))); ...
                    sum(sum(Hessmix(1:self.sng1,:))) sum(sum(Hess2(1:self.sng1,1:self.sng1))) sum(sum(Hess2(self.sng1+1:end,1:self.sng1))); ...
                    sum(sum(Hessmix(self.sng1+1:end,:))) sum(sum(Hess2(self.sng1+1:end,1:self.sng1))) sum(sum(Hess2(self.sng1+1:end,self.sng1+1:end)))];
        end
        function Plot(self,hest)
            %This can be used to plot the fitting of the data during the
            %solving of the homogeneous estimate.
            self.solver.sigma = [hest(2)*ones(self.sng1,1); hest(3)*ones(self.ng-self.sng1,1)];
            self.solver.Plot(hest(1)*ones(self.ngc,1));
        end
        
        
    end
    
end
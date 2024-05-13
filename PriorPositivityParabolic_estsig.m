classdef PriorPositivityParabolic_estsig < handle
    %PositivityParabolic class contains a simple parabolic positivity
    %constraint and the required functions to be called by an inverse
    %problems solver (e.g. SolverLinesearch.m).
    %Supports real values, complex values, or a stacked complex model
    %(where input is [real(sigma); imag(sigma)])
    %
    %Author: Petri Kuusela 16.11.2021
    properties
        a %the constraint function value will be a*(x-minAcceptable)^2 for all x < minAcceptable
        minAcceptable%The minimum value of estimate that is not penalized
    end
    methods
        function obj = PriorPositivityParabolic_estsig(minAcceptable, valAt0)
            %Class constructor.
            %Input: minAcceptable = The minimum value of estimate which is
            %                       not penalized.
            %       valAt0 = What is the functional value at 0
            if isstruct(minAcceptable)
                obj.minAcceptable = minAcceptable;
            else
                obj.minAcceptable.sigma = minAcceptable;
                obj.minAcceptable.coupl = minAcceptable;
                obj.minAcceptable.zeta = minAcceptable;
            end                                 
            if isstruct(valAt0)              
                obj.a.sigma = valAt0.sigma./obj.minAcceptable.sigma;
                obj.a.coupl = valAt0.coupl./obj.minAcceptable.coupl;
                obj.a.zeta = valAt0.zeta./obj.minAcceptable.zeta;
            else
                obj.a.sigma = valAt0./obj.minAcceptable.sigma;
                obj.a.coupl = valAt0./obj.minAcceptable.coupl;
                obj.a.zeta = valAt0./obj.minAcceptable.zeta;
            end
        end
        function res = OptimizationFunction(self, est)
            res = 0;
            if ~isempty(est.sigma)
                res = self.a.sigma*sum((est.sigma(est.sigma<self.minAcceptable.sigma) - self.minAcceptable.sigma).^2);
            end
            if ~isempty(est.coupl)
                res = res + self.a.coupl*sum((est.coupl(est.coupl<self.minAcceptable.coupl) - self.minAcceptable.coupl).^2);
            end
            if ~isempty(est.zeta)
                res = res + self.a.zeta*sum((est.zeta(est.zeta<self.minAcceptable.zeta) - self.minAcceptable.zeta).^2);
            end
        end
        function [Hess, grad] = GetHessAndGrad(self, est)            
            if ~isempty(est.sigma)
                grad.sigma = zeros(length(est.sigma),1);
                select = est.sigma<self.minAcceptable.sigma;
                grad.sigma(select) = 2*self.a.sigma*(est.sigma(select) - self.minAcceptable.sigma);
                Hessvec = zeros(length(est.sigma),1);
                Hessvec(select) = 2*self.a.sigma;%./self.minAcceptable^2;
                Hess.sigma = diag(Hessvec);
            else
                grad.sigma = [];
                Hess.sigma = [];
            end
            if ~isempty(est.coupl)
                grad.coupl = zeros(length(est.coupl),1);
                select = est.coupl<self.minAcceptable.coupl;
                grad.coupl(select) = 2*self.a.coupl*(est.coupl(select) - self.minAcceptable.coupl);
                Hessvec = zeros(length(est.coupl),1);
                Hessvec(select) = 2*self.a.coupl;%./self.minAcceptable^2;
                Hess.coupl = diag(Hessvec);
            else
                grad.coupl = [];
                Hess.coupl = [];
            end
            if ~isempty(est.zeta)
                grad.zeta = zeros(length(est.zeta),1);
                select = est.zeta<self.minAcceptable.zeta;
                grad.zeta(select) = 2*self.a.zeta*(est.zeta(select) - self.minAcceptable.zeta);
                Hessvec = zeros(length(est.zeta),1);
                Hessvec(select) = 2*self.a.zeta;%./self.minAcceptable^2;
                Hess.zeta = diag(Hessvec);
            else
                grad.zeta = [];
                Hess.zeta = [];
            end

            grad = DLSS_est(grad.sigma, grad.coupl, grad.zeta);
            Hess = DLSS_Hess(Hess.sigma, Hess.coupl, Hess.zeta);
        end
        
    end
end
        
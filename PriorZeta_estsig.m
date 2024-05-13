classdef PriorZeta_estsig < handle
    %Author: Petri Kuusela 19.8.2023
    properties
        ev  %expected value
        a   %strength
    end
    methods
        function obj = PriorZeta_estsig(ev, a)
            obj.ev = ev;
            obj.a = a;
        end
        function res = OptimizationFunction(self, est)
            res = 0;
            if ~isempty(est.zeta)
                res = res + self.a*sum((est.zeta - self.ev).^2);
            end
        end
        function [Hess, grad] = GetHessAndGrad(self, est)            
            if ~isempty(est.zeta)
                grad = 2*self.a*(est.zeta - self.ev);
                Hessvec = 2*self.a*ones(length(est.zeta),1);
                Hess = diag(Hessvec);
            else
                grad = [];
                Hess = [];
            end
            if ~isempty(est.sigma)
                sgrad = zeros(size(est.sigma));
                sHess = zeros(size(est.sigma,1));
            else
                sgrad = [];
                sHess = [];
            end
            if ~isempty(est.coupl)
                cgrad = zeros(size(est.coupl));
                cHess = zeros(size(est.coupl,1));
            else
                cgrad = [];
                cHess = [];
            end

            grad = DLSS_est(sgrad, cgrad, grad);
            Hess = DLSS_Hess(sHess, cHess, Hess);
        end
        
    end
end
        
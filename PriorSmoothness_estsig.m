classdef PriorSmoothness_estsig < handle
%PriorSmoothness class contains a smoothness prior and the required functions
%for it to be called by an inverse problem solver -class (e.g.
%SolverLinesearch.m).
%
%Supports real values, complex values, or a stacked complex model
%(where input is [real(sigma); imag(sigma)]). To enable using stacked
%model, the mean has to have 2 elements or twice the length of ginv.
%
%Supports 2D and 3D meshes (3D has not been tested though...).
%
%Author: Petri Kuusela 16.11.2021
    
    
    properties
        invcov  %The inverse of the covariance matrix
        mean    %The expectation value of the estimate (can be a single value or a vector with same length as ginv)
        L       %chol(invcov), used only to draw samples from the prior
        regularization_c %A regularization constant to make sure the covariance matrix is not poorly conditioned
        common_c %A constant added to all the elements of the covariance matrix. Larger values allow greater variation in the mean value of the estimate
    end
    
    
    methods
        
        function obj = PriorSmoothness_estsig(ginv, corlen, var, mean)
            %Class constructor.
            %Input: ginv = the nodes of the reconstruction mesh
            %       corlen = the correlation length (length where the
            %               cross-covariance has decreased to 1%)
            %       var = the variance of the values of the estimate
            %       mean = the expectation value of the estimate. Can be a
            %       single value or a vector with same length as ginv for a
            %       non-stacked model, and two values or twice the length
            %       of ginv for a stacked model (i.e. where the estimate is
            %       [real(sigma);imag(sigma)])
            
            obj.regularization_c = 1e-4;
            obj.common_c = 1e-1;
            obj.SetCovMat(ginv, corlen, var);
            if isstruct(mean)
                obj.mean = mean;
            else
                obj.mean.sigma = mean;
                obj.mean.coupl = mean;
            end
            obj.L.sigma = chol(obj.invcov.sigma);
            obj.L.coupl = chol(obj.invcov.coupl);

        end
        
        function SetCovMat(self, g, corlen, var)
            %Compute the inverse covariance matrix given the input
            %parameters. 
            %Input: g = node coordinates
            %       corlen = spatial correlation length (length where the
            %               cross-covariance has decreased to 1%)
            %       var = the variance of the estimate values
            %       stacked = A flag: create a stacked model?
            
            if ~isstruct(g)
                tg = g;
                clear g;
                g.sigma = tg;
                g.coupl = tg;
            end
            if ~isstruct(corlen)
                tcorlen = corlen;
                clear corlen;
                corlen.sigma = tcorlen;
                corlen.coupl = tcorlen;
            end
            if ~isstruct(var)
                tvar = var;
                clear var;
                var.sigma = tvar;
                var.coupl = tvar;
            end
            
            ng = size(g.sigma,1);
            b = corlen.sigma./sqrt(2*log(100));
            c = self.regularization_c*var.sigma;
            a = var.sigma-c;
            xmat = repmat(g.sigma(:,1),1,ng);
            ymat = repmat(g.sigma(:,2),1,ng);
            zmat = 0;
            if size(g.sigma,2) > 2%if the mesh is 3D
                zmat = repmat(g.sigma(:,2),1,ng);
            end
            cov = a*exp(-0.5*((xmat-xmat').^2+(ymat-ymat').^2+(zmat-zmat').^2)/b^2);
            cov = cov + diag(c*ones(size(cov,1),1)) + self.common_c*var.sigma*ones(size(cov,1));
            self.invcov.sigma = inv(cov);

            ng = size(g.coupl,1);
            b = corlen.coupl./sqrt(2*log(100));
            c = self.regularization_c*var.coupl;
            a = var.coupl-c;
            xmat = repmat(g.coupl(:,1),1,ng);
            ymat = repmat(g.coupl(:,2),1,ng);
            zmat = 0;
            if size(g.coupl,2) > 2%if the mesh is 3D
                zmat = repmat(g.coupl(:,2),1,ng);
            end
            cov = a*exp(-0.5*((xmat-xmat').^2+(ymat-ymat').^2+(zmat-zmat').^2)/b^2);
            cov = cov + diag(c*ones(size(cov,1),1)) + self.common_c*var.coupl*ones(size(cov,1));
            self.invcov.coupl = inv(cov);
        end
        
        function res = OptimizationFunction(self, est)
            %This function is called by the inverse problem solver, and it
            %gives the function value to be minimized.
            res = 0;
            if ~isempty(est.sigma)
                res = res + 0.5*(est.sigma-self.mean.sigma)'*self.invcov.sigma*(est.sigma-self.mean.sigma);
            end
            if ~isempty(est.coupl)
                res = res + 0.5*(est.coupl-self.mean.coupl)'*self.invcov.coupl*(est.coupl-self.mean.coupl);
            end
        end
        
        function [Hess, grad] = GetHessAndGrad(self, est)
            %This function is called by the inverse problem solver, and it
            %gives the Hess-matrix and gradient of the optimization
            %function.
            if ~isempty(est.sigma)
                grads = self.invcov.sigma*(est.sigma-self.mean.sigma);
                Hesss = self.invcov.sigma;
                if size(est.sigma,1) == 1
                    grads = sum(grads);
                    Hesss = sum(sum(Hesss));
                end
            else
                grads = [];
                Hesss = [];
            end
            if ~isempty(est.coupl)
                gradc = self.invcov.coupl*(est.coupl-self.mean.coupl);
                Hessc = self.invcov.coupl;
                if size(est.coupl,1) == 1
                    gradc = sum(gradc);
                    Hessc = sum(sum(Hessc));
                end
            else
                gradc = [];
                Hessc = [];
            end
            grad = DLSS_est(grads, gradc, zeros(size(est.zeta)));
            Hess = DLSS_Hess(Hesss, Hessc, zeros(size(est.zeta,1)));
        end


        function samples = DrawSamples(self, n)
            vec = randn(size(self.L.sigma,1),n);
            if length(self.mean.sigma) == 1
                samples.sigma = self.L.sigma\vec + self.mean.sigma*ones(size(self.L.sigma,1),n);
            else
                samples.sigma = self.L.sigma\vec + self.mean.sigma*ones(1,n);
            end

            vec = randn(size(self.L.coupl,1),n);
            if length(self.mean.coupl) == 1
                samples.coupl = self.L.coupl\vec + self.mean.coupl*ones(size(self.L.coupl,1),n);
            else
                samples.coupl = self.L.coupl\vec + self.mean.coupl*ones(1,n);
            end

        end


        function d = ComputeCuttingLength(~, ginv, n)
            %This is currently not in use and has not been tested
            sn = n/size(ginv,1);
            avgsize = 1/size(ginv,2)*(sum(max(ginv)) - sum(min(ginv)));
            d = (sn/size(ginv,1))^(1/size(ginv,2))*avgsize;
        end
        
        
    end
    
    
    
end
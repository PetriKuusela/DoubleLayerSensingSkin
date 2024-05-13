classdef PriorTotalVariation_estsig
%PriorTotalVariation class contains a TV prior and the required functions
%for it to be called by an inverse problem solver -class (e.g.
%SolverLinesearch.m).
%Supports real values, complex values, or a stacked complex model
%(where input is [real(sigma); imag(sigma)])
%
%Author: Petri Kuusela 16.11.2021
    
    properties
        gradphis     %The gradients of all basis functions. A cell array with each cell containing gradients inside one element
        gradphic
        Hs           %The elements of the reconstruction mesh
        gs           %The nodes of the reconstruction mesh
        Hc
        gc
        nHs          %Number of elements
        ngs          %Number of nodes
        nHc
        ngc
        alpha       %The "strength" coefficient of TV
        beta        %The smoothing parameter
        cf          %Coupling factor, how strongly real and imag parts couple
        Areass       %A vector containing the area of each element
        Areasc       %A vector containing the area of each element
        D           %This is required only for Chambolle-Pock type solvers
        gdim        %dimension of the mesh
    end
    
    methods
        
        function obj = PriorTotalVariation_estsig(gs, Hs, gc, Hc, ec, ef, beta)
            %Class constructor.
            %Inputs: g = the nodes of the reconstruction mesh
            %        H = the elements of the reconstruction mesh(for these
            %        2: s = sigma, c = coupl)
            %        ec = expected change of the estimate (can have two
            %        elements if different coefficients are to be used for
            %        conductivity and permittivity.)
            %        ef(optional) = extra factor which multiplies alpha
            %        beta(optional) = the smoothing constant
            %NOTE: if you want to set alpha yourself, you can change it
            %directly after initiating the object.
            %NOTE2: Older version of this class constructor took in alpha
            %as an argument, so that usage may be found in older drivers.
            
            if nargin < 6 || isempty(ef)
                ef = 1;
            end
            if nargin < 7 || isempty(beta)
                beta = 1e-4;
            end

            if ~isstruct(ec)
                tec = ec;
                clear ec;
                ec.sigma = tec;
                ec.coupl = tec;
            end
            if ~isstruct(ef)
                tef = ef;
                clear ef;
                ef.sigma = tef;
                ef.coupl = tef;
            end
            if ~isstruct(beta)
                teb = beta;
                clear beta;
                beta.sigma = teb;
                beta.coupl = teb;
            end
            
            obj.gs = gs;
            obj.Hs = Hs;
            obj.gc = gc;
            obj.Hc = Hc;

            mesh_a = prod(max(gs)-min(gs))/size(Hs,1);%roughly the element area/volume of the mesh
            mesh_c = (2*mesh_a)^(1/size(gs,2));%roughly the element size of mesh
            egrad = ec.sigma./mesh_c;%expected gradient
            tvp = 0.975;
            alpha = -ef.sigma.*log(1-tvp)./(mesh_a*egrad);
            obj.alpha.sigma = alpha;

            mesh_a = prod(max(gc)-min(gc))/size(Hc,1);%roughly the element area/volume of the mesh
            mesh_c = (2*mesh_a)^(1/size(gc,2));%roughly the element size of mesh
            egrad = ec.coupl./mesh_c;%expected gradient
            tvp = 0.975;
            alpha = -ef.coupl.*log(1-tvp)./(mesh_a*egrad);
            obj.alpha.coupl = alpha;

            obj.beta = beta;
            
            obj.gdim = size(gs,2);
            obj.ngs = size(gs, 1);
            obj.nHs = size(Hs, 1);
            obj.ngc = size(gc, 1);
            obj.nHc = size(Hc, 1);
            [obj.gradphis, obj.Areass, obj.D] = obj.ComputeGrads2D(gs,Hs);
            [obj.gradphic, obj.Areasc, obj.D] = obj.ComputeGrads2D(gc,Hc);
            
        end
        
        
        function res = OptimizationFunction(self, est)
            %COMMENTS OUTDATED!!!
            %This is the function called by the inverse problem solver to
            %get the value of the function which is to be minimized.
            %Input: sigest is the conductivity, can be either real or
            %       complex, or stacked so that real and imaginary part are
            %       on top of each other.
        
            gradsigma = self.ComputeGrads(est);
            res = 0;
            if ~isempty(est.sigma)
                res = self.alpha.sigma*sum(self.Areass.*sqrt(sum(abs(gradsigma.sigma).^2, 2)+self.beta.sigma));
            end
            if ~isempty(est.coupl)
                res = res + self.alpha.coupl*sum(self.Areasc.*sqrt(sum(abs(gradsigma.coupl).^2, 2)+self.beta.coupl));
            end
        end
        
        function [Hess, grad] = GetHessAndGrad(self, est)
            %COMMENTS OUTDATED, AS WITH ALL ESTSIGs!!!
            %Compute the Hess-matrix and gradient of the optimization
            %function. This function is called by the inverse problem
            %solver.
            %Input: sigest is the conductivity, can be either real or
            %       complex, or stacked so that real and imaginary part are
            %       on top of each other.
            
            a = self.alpha;
            gradest = self.ComputeGrads(est);%gradsigma is a cell-array to facilitate similart treatment as in complex case
            
            grads = [];
            gradc = [];
            Hesss = [];
            Hessc = [];

            if ~isempty(est.sigma)
                grads = zeros(self.ngs,1);%initialize arrays for gradient and Hess-matrix
                Hvals = zeros(6*self.ngs,1);%Hess matrix is collected in sparse form
                HindsI = zeros(6*self.ngs,1);
                HindsJ = zeros(6*self.ngs,1);
                
                II = 1;%index for collecting the Hess-values in sparce matrix       
                %This block calculates the gradient of the TV functional
                for ii=1:self.nHs
                    tdgrads = self.gradphis{ii};%The gradients of the basis functions in element ii
                    tinds = self.Hs(ii,:);%node indices that are the verices of element ii
                    tdsigma = gradest.sigma(ii,:);%The gradient of sigma in element ii
                    gnorm = 1/sqrt(sum(tdsigma.^2) + self.beta.sigma);%1/sqrt part of the gradient
                    tdot = tdgrads*tdsigma';%The dot product of gradients of basis functions and gradient of sigma
                    grads(tinds) = grads(tinds) + a.sigma*gnorm.*tdot.*self.Areass(ii);%Add the terms to the relevant gradient elements
                end
    
                %This block calculates the Hess-matrix of the TV functional
                for ii=1:self.nHs
                    tdgrads = self.gradphis{ii};%The gradients of the basis functions in element ii
                    tinds = self.Hs(ii,:);%node indices that are the vertices of element ii
                    tdsigma = gradest.sigma(ii,:);%The gradient of sigma in element ii
                    gnorm = (sum(tdsigma.^2) + self.beta.sigma)^(-0.5);%These are parts of the Hessian
                    gnorm2 = (sum(tdsigma.^2) + self.beta.sigma)^(-1.5);
                    
                    for jj=1:3
                        for kk=1:3
                            phii = tdgrads(jj,:);%grad? of basis function_i
                            phij = tdgrads(kk,:);%grad? of basis function_j
                            f = phii*phij'*gnorm;
                            h = -(tdsigma*phii'*gnorm2*tdsigma*phij');
                            Hvals(II) = a.sigma*(f + h).*self.Areass(ii);
                            HindsI(II) = tinds(jj);
                            HindsJ(II) = tinds(kk);
                            II = II + 1;
                        end
                    end
                end
                Hesss = sparse(HindsI,HindsJ, Hvals);

            end

            %BEGIN COUPLING PART:

            if ~isempty(est.coupl)

                gradc = zeros(self.ngc,1);%initialize arrays for gradient and Hess-matrix
                Hvals = zeros(6*self.ngc,1);%Hess matrix is collected in sparse form
                HindsI = zeros(6*self.ngc,1);
                HindsJ = zeros(6*self.ngc,1);
                
                II = 1;%index for collecting the Hess-values in sparce matrix       
                %This block calculates the gradient of the TV functional
                for ii=1:self.nHc
                    tdgrads = self.gradphic{ii};%The gradients of the basis functions in element ii
                    tinds = self.Hc(ii,:);%node indices that are the verices of element ii
                    tdsigma = gradest.coupl(ii,:);%The gradient of sigma in element ii
                    gnorm = 1/sqrt(sum(tdsigma.^2) + self.beta.coupl);%1/sqrt part of the gradient
                    tdot = tdgrads*tdsigma';%The dot product of gradients of basis functions and gradient of sigma
                    gradc(tinds) = gradc(tinds) + a.coupl*gnorm.*tdot.*self.Areasc(ii);%Add the terms to the relevant gradient elements
                end
    
                %This block calculates the Hess-matrix of the TV functional
                for ii=1:self.nHc
                    tdgrads = self.gradphic{ii};%The gradients of the basis functions in element ii
                    tinds = self.Hc(ii,:);%node indices that are the vertices of element ii
                    tdsigma = gradest.coupl(ii,:);%The gradient of sigma in element ii
                    gnorm = (sum(tdsigma.^2) + self.beta.coupl)^(-0.5);%These are parts of the Hessian
                    gnorm2 = (sum(tdsigma.^2) + self.beta.coupl)^(-1.5);
                    
                    for jj=1:3
                        for kk=1:3
                            phii = tdgrads(jj,:);%grad? of basis function_i
                            phij = tdgrads(kk,:);%grad? of basis function_j
                            f = phii*phij'*gnorm;
                            h = -(tdsigma*phii'*gnorm2*tdsigma*phij');
                            Hvals(II) = a.coupl*(f + h).*self.Areasc(ii);
                            HindsI(II) = tinds(jj);
                            HindsJ(II) = tinds(kk);
                            II = II + 1;
                        end
                    end
                end
                Hessc = sparse(HindsI,HindsJ, Hvals);

            end

            %grad = DLSS_est(grads, gradc, zeros(size(est.zeta)));
            %Hess = DLSS_Hess(Hesss, Hessc, zeros(size(est.zeta,1)));
            grad = DLSS_est(grads, gradc, zeros(size(est.zeta)));
            Hess = DLSS_Hess(Hesss, Hessc, zeros(size(est.zeta,1)));

            
        end
        
        function gradest = ComputeGrads(self, est)
            %Compute the gradients of the sigma, using the pre-computed
            %self.gradphi.
            gradest = [];
            if ~isempty(est.sigma)
                N = size(self.gradphis,1);
                gradest.sigma = zeros(N,self.gdim);
                for ii=1:N
                    sigmas = est.sigma(self.Hs(ii,:));
                    grads = self.gradphis{ii};
                    gradest.sigma(ii,:) = sum(diag(sigmas)*grads); %sigmas(1)*grads(1,:) +  sigmas(2)*grads(2,:) +  sigmas(3)*grads(3,:);
                end
            end

            if ~isempty(est.coupl)
                N = size(self.gradphic,1);
                gradest.coupl = zeros(N,self.gdim);
                for ii=1:N
                    coupls = est.coupl(self.Hc(ii,:));
                    grads = self.gradphic{ii};
                    gradest.coupl(ii,:) = sum(diag(coupls)*grads); %sigmas(1)*grads(1,:) +  sigmas(2)*grads(2,:) +  sigmas(3)*grads(3,:);
                end
            end

        end
        
        function [gradphi, areas, R] = ComputeGrads2D(self, g,H)
            % Computes the spatial gradients of the linear basis functions
            % Also calculates the area of each element, and R, which is
            % only needed for Chambolle-Pock type solvers.
            N = size(H,1);
            gN = size(g,1);
            gradphi = cell(N,1);
            areas = zeros(N,1);
            gradsigma = zeros(N,2);
            L = [-ones(size(g,2),1) eye(size(g,2))];
            R = sparse(2*N,gN);
            for ii=1:N
                tind = H(ii,:);
                X = g(H(ii,:),:);
                Jt = L*X;
                areas(ii) = 1/2*abs(det(Jt));
                grads = Jt\L;
                grads = grads';
                gradphi{ii} = grads;
                %R(ii,tind) = grads(:,1)';%These are only needed for Chambolle-Pock type solvers
                %R(end/2 + ii,tind) = grads(:,2)';%Should these be multiplied by areas? 
                R = 0;%This is not working currently
            end

        end


        %%
        %The rest are for Chambolle-Pock type inverse problem solvers and
        %are most likely outdated
        
        function D = matrix(self, sigma)
            D = self.D;
        end
        function res = Fproximal(self, d, s)
            dl = sqrt(d(1:end/2).^2 + d(end/2+1:end).^2);
            dl = [dl;dl];
            res = self.alpha*self.Areasx2.*d./max(self.alpha*self.Areasx2, dl);      
        end
        function res = nextDual(self, d, s, sigest2)
            res = d + self.matrix(sigest2)*sigest2*s;
        end 
        
        
    end 
    
end
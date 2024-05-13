classdef EITFEM_DLSS < EITFEM
    
    
    properties
        sigma
        Jr
        Aic
        Avc
        couplforJ
        mincoupl
    end
        
    methods
        
        function obj = EITFEM_DLSS(fmesh)
            
            obj@EITFEM(fmesh);
            obj.Jr = fmesh.RhsJacobianTerms();
            [obj.Aic, obj.Avc] = fmesh.GradientMatrixCoupling();
            obj.mincoupl = 1e-6;
            
        end
        
        function elval = SolveForward(self, sigma, coupling, noscf)
            %Solve the potentials or currents for given sigma and angular
            %frequency omega.
            %Input:     sigma   = The conductivity distribution
            %           mode    = 'current' or 'potential', which one is
            %                     injected
            %                     
            %           noscf   = Do not use scaling of sigma
            %Output: elval =  the currents or potentials for a single
            %                 frequency in matrix (and complex) form
            %                 computed using FEM  
            
            
            if nargin < 4 || isempty(noscf)
                noscf = 0;
            end
            if length(coupling) == 1 && coupling == 1%this is a bad hack due to bad design. coupling == 1 means this is called from Jacobian subroutine
                noscf = 1;
                coupling = self.couplforJ;
            end
            
                
            %Check formats of sigma and self.scales and lift any
            %values that are negative or too close to 0.
            sigma = self.LiftSigma(sigma);
            coupling(coupling<self.mincoupl) = self.mincoupl;
            
            %Scale sigma to avoid numerical problems:
            if noscf
                scf = 1;%Use scf=1 for updating the intermidiate steps  used in calculating Jacobian
            else
                scf = 1;%Use of scf disabled due to it conctantly causing bugs
                %scf = abs(mean(sigma));%This is used to avoid some numerical problems following from poorly scaled sigma.
            end
            sigma = sigma./scf;%This is not the optimal scaling, and may not work in some extreme cases
            zetas = self.zeta.*scf;
            coupling = coupling./scf;
            
            
            %Start forming the EIT-matrix
            [A0, rhs] = self.fmesh.SigmadPhiidPhij(sigma, coupling, reshape(self.IDmeas, self.fmesh.nel, numel(self.IDmeas)/self.fmesh.nel));
            if strcmp(self.mode, 'current')
                error('current injection is not yet implemented');
            end
            
            if self.stackOutput
                Inj = self.IDmeas(1:end/2) + 1i*self.IDmeas(end/2+1:end);
            else
                Inj = self.IDmeas;
            end
            if self.recalc
                self.S = self.intS{1}/zetas(1);
                for ii=2:length(self.intS)
                    self.S = self.S + self.intS{ii}/zetas(ii);
                end
                if strcmp(self.mode, 'potential')
                    tempM = self.intM*diag(1./zetas);
                    self.S = [self.S zeros(size(self.S,1), size(self.C,2));...
                              -self.C'*tempM' self.C'*self.C];
                    ninj = numel(Inj)/(self.fmesh.nel);
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(Inj, self.fmesh.nel, ninj);
                    self.b(1:self.fmesh.ng,:) = tempM*InjMat;
                    self.b(end-self.fmesh.nel+2:end,:) = ...%self.b(end-self.fmesh.nel+2:end,:)...
                                        - self.C'*(diag(self.intB./zetas)*InjMat);
                elseif strcmp(self.mode, 'current')%current injection is not yet implemented!
                    C2 = -self.C'*diag(1./zetas)*self.intM';
                    self.S = [self.S C2'; C2 -self.C'*diag(self.intB./zetas)*self.C];
                    ninj = numel(Inj)/self.fmesh.nel;
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(Inj, self.fmesh.nel, ninj);
                    self.b(end-ninj+2:end,:) = self.C'*InjMat;
                else
                    error(['Unrecognized solver mode: ' self.mode]);
                end
                self.recalc = 0;
            end
            
            self.A = A0 + self.S;

            self.Pot = self.A\(self.b + rhs);

            elval = self.QC*self.Pot;
            if strcmp(self.mode, 'potential')
                elval = elval*scf;%The results have to be scaled in the opposite direction as the conductivity was scaled in the beginning
            elseif strcmp(self.mode, 'current')
                elval = elval/scf;%The results have to be scaled in the same direction as the conductivity was scaled in the beginning
            else
                error(['Unrecognized solver mode: ' self.mode]);
            end
            
        end %end solveForward
        
        function vec = SolveForwardVec(self, sigma, coupling)
            %Solve FEM with given sigma.
            %output: vec = the currents or potentials in vector format
            %              computed using FEM
            
            if nargin == 2 %assume "sigma" is actually "coupling"
                coupling = sigma;
                sigma = self.sigma;
            end
            
            vec = self.SolveForward(sigma, coupling);
            vec = vec(:);
            vec = vec(self.vincl);
            vec = vec + self.eps;%Add a given constant (scalar or vector) to the results
            if self.stackOutput
                vec = [real(vec); imag(vec)];
            end
        end %end solveForwardVec
        
        function res = OptimizationFunction(self, coupl)
            %Calculate the residual of the forward problem, which is to be
            %minimized (in addition to the regularization) when solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)        
           
            elval = self.SolveForwardVec(self.sigma, coupl);
            res = 0.5*sum((self.Ln*(abs(self.Dmeas - elval))).^2);
        end
        
        function [Hess, grad, Hessmix] = GetHessAndGradSigma(self, coupl, sigma)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)
            
            
            self.sigma = sigma;
            J = self.FullJacobian(sigma);
            Hess = J'*self.InvGamma_n*J;%First order approximation of the forward problem
            fres = self.SolveForwardVec(sigma, coupl);
            grad = J'*self.InvGamma_n*(fres - self.Dmeas);
            Jc = self.FullJacobianCoupl(coupl);
            Hessmix = J'*self.InvGamma_n*Jc;
            
        end
        
        function [Hess, grad] = GetHessAndGrad(self, coupling)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)
            
            
            J = self.FullJacobianCoupl(coupling);
            Hess = J'*self.InvGamma_n*J;%First order approximation of the forward problem
            fres = self.SolveForwardVec(coupling);
            grad = J'*self.InvGamma_n*(fres - self.Dmeas);
            
        end
        
        function J = FullJacobianCoupl(self, coupl)
            %Calculate the Jacobian of the forward problem at sigma            
            
            J = self.JacobianCoupl(coupl);
            if length(self.scales) == 2%The scaling done before solving the FEM has to be taken into account here.
                J(:,1:end/2) = J(:,1:end/2)*self.scales(1);
                J(:,end/2+1:end) = J(:,end/2+1:end)*self.scales(2);
            else
                J = J*diag(self.scales);%The scaling done before solving the FEM has to be taken into account here.
            end
            
        end
        
        function J = JacobianCoupl(self, coupl)
            %Calculate the Jacobian of the forward problem at sigma with
            %injection frequency omega. For info on mode and meshtype see
            %function SolveForward
            
            elval = self.SolveForward(self.sigma, coupl, 1);%This elval may be inaccurate since noscf=1, but the values are not used anywhere
            %if exist('Jacobian3dECTmx.mexw64')
            if 0==1%Using the mex function is currently disabled for easier debugging
                Jleft = -full(self.QC/self.A);%Using the mex function requires some pre-processing
                tempStruct.C = self.QC;
                tempStruct.Electrode = Iel;
                tempStruct.Pot = self.Pot(1:self.ng, :);
                tempStruct.PotF = self.Pot;
                Js = Jacobian3dECTmx(self.Ai, self.Av, Jleft, tempStruct);
            else
                Js = self.Jacobian3dECTCoupl(elval);
            end
            Js = Js(self.vincl,:);%some post processing of the Jacobian
            J = self.fmesh.JacobianFtoIc(Js);
            
            
            %The final output format of J is governed mainly by three
            %variables in self: absMeas, stackOutput, and self.nginv (or
            %length(sigma) compared to that, i.e. is input stacked?).
            
            if self.stackOutput
                if length(coupl) == 2*self.fmesh.nginv%a complex stacked model
                    J = [real(J) -imag(J); imag(J) real(J)];
                elseif any(~isreal(coupl))% a complex non-stacked model (TODO this is not yet implemented)
                    error('Non-stacked complex sigma is not supported with stacked output.\nThis would require the stacking done after calculating the gradient.');
                else %only real values in sigma, but still want stacked model?
                    warning('real valued sigma, but complex stacked output required, are you sure?');
                    J = [J; zeros(size(J))];
                end
            end
        end
                
        function Js = Jacobian3dECTCoupl(self, elval)
            % Computes the Jacobian J = d(measurements) / d(coupl)
            % Original CEM-ECT version written by K. Karhunen 03.06.2013 is implementation
            % from the earlier versions by V. Rimpiläinen(?) and L. M. Heikkinen.
            % Current version adapted for use by P. Kuusela 15.11.2021
            b = length(self.Aic);
            c = size(elval,1); 
            d = size(self.Pot,2); 
            
            Uel = reshape(self.IDmeas, c, d);
            Jleft  = self.QC/self.A;

            Jright = -self.Pot;

            Js = zeros(c*d,b);

            for ii=1:b
              Jid = self.Aic{ii};

              Jtemp   = Jleft(:,Jid)*self.Avc{ii}*Jright(Jid,:);
              Js(:,ii) = Jtemp(:);
            end
            for ii = 1:length(self.Jr)
                if ~isempty(self.Jr{ii})
                    Jtemp = Jleft*self.Jr{ii}*Uel(self.fmesh.elng(ii),:);
                    Js(:,ii) = Js(:,ii) + Jtemp(:);
                end
            end

        end
                
        function SetInvGammaCorrelated(self, basenoise, elnoise)
            %Old Comments, not up to date!!!
            %Calculate and set the inverse of covariance matrix based on
            %the noise levels given as arguments.
            %Input: meas_noise_coef_e = a constant noise level for all
            %       measurements. This coefficient is scaled to the
            %       difference of minimum and maximum measurement values to
            %       get the actual noise level. Can be a scalar or vector
            %       of length self.omegas having a different coefficient
            %       for different frequencies. A complex value can be used
            %       to determine different values for real and imag parts
            %       of measurements.
            %
            %       meas_noise_coef2(optional) = relative noise level of the
            %       measurements. This is multiplied by the absolute value
            %       of each measurement to get the noise level of that
            %       measurement. Can be a scalar or vector
            %       of length self.omegas having a different coefficient
            %       for different frequencies. A complex value can be used
            %       to determine different values for real and imag parts
            %       of measurements.
            %
            %       The two types of noise above are added together.
            %NOTE: This does not add noise to the solution of the forward
            %problem! This just calculates the Weighing matrix used in the
            %inverse problem.
            
            nmeas = numel(self.Dmeas);            
            
            %same coefficients for real and imag parts
            %select = self.vincl;%1:nmeas;
            %Calculating the model variance of the noise
            %Constant noise for all measurements:
            var_meas = (basenoise*(max(max(abs(self.Dmeas)))-min(min(abs(self.Dmeas)))))^2;
            %Noise level relative to the abs of measurement value:
            %var_meas = var_meas + (meas_noise_coef2*(abs(self.Dmeas) - min(min(abs(self.Dmeas))))).^2;
            %Assume no cross-correlations, hence a diagonal covariance mat:
            Gamma_n = diag(var_meas(:));
            self.InvGamma_n = sparse(inv(Gamma_n));
            self.Ln = chol(self.InvGamma_n);%Store both invGamma and it's Cholesky

        end
        
        
    end
    
    
    
    
end
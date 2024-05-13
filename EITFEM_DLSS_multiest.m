classdef EITFEM_DLSS_multiest < EITFEM


    properties
        sigma
        Jr
        Aic
        Avc
        coupl
        mincoupl
        minzeta
    end
        
    methods
        
        function obj = EITFEM_DLSS_multiest(fmesh)
            
            obj@EITFEM(fmesh);
            obj.Jr = fmesh.RhsJacobianTerms();
            [obj.Aic, obj.Avc] = fmesh.GradientMatrixCoupling();
            obj.mincoupl = 1e-6;
            obj.minzeta = 1e-9;
            
        end
        
        function elval = SolveForward(self, est, noscf)
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
            
            %WARNING: Horrible hack ahead!
            if ~isstruct(est) && ~isobject(est)
                tempvar = est;
                clear est;
                est.sigma = tempvar;
                est.coupl = [];
                est.zeta = [];
            end
            
            if nargin < 3 || isempty(noscf)
                noscf = 1;
            end
            if ~isempty(est.coupl)
                if size(est.coupl,1) > 1
                    tcoupl = est.coupl;
                else
                    tcoupl = est.coupl*ones(size(self.fmesh.itofc,2),1);
                end
            else
                tcoupl = self.coupl;
            end
            if ~isempty(est.sigma)
                if size(est.sigma,1) > 1
                    tsigma = est.sigma;
                else
                    tsigma = est.sigma*ones(size(self.fmesh.itof,2),1);
                end
            else
                tsigma = self.sigma;
            end
            if ~isempty(est.zeta)
                if size(est.zeta,1) > 1
                    tzeta = est.zeta;
                else
                    tzeta = est.zeta*ones(self.fmesh.nel,1);
                end
            else
                tzeta = self.zeta;
            end
            %if length(coupling) == 1 && coupling == 1%this is a bad hack due to bad design. coupling == 1 means this is called from Jacobian subroutine
            %    noscf = 1;
            %    coupling = self.couplforJ;
            %end
            
            %Check formats of sigma and self.scales and lift any
            %values that are negative or too close to 0.
            tsigma = self.LiftSigma(tsigma);
            tcoupl(tcoupl<self.mincoupl) = self.mincoupl;
            tzeta(tzeta<self.minzeta) = self.minzeta;
            
            %Scale sigma to avoid numerical problems:
            if noscf
                scf = 1;%Use scf=1 for updating the intermidiate steps  used in calculating Jacobian
            else
                %Use of scf is discouraged for it continuously causing bugs
                scf = abs(mean(tsigma));%This is used to avoid some numerical problems following from poorly scaled sigma.
            end
            tsigma = tsigma./scf;%This is not the optimal scaling, and may not work in some extreme cases
            tzeta = tzeta.*scf;
            tcoupl = tcoupl./scf;
            
            
            %Start forming the EIT-matrix
            [A0, rhs, Cc, elcs] = self.fmesh.SigmadPhiidPhij(tsigma, tcoupl, reshape(self.IDmeas, self.fmesh.nel, numel(self.IDmeas)/self.fmesh.nel));
            if strcmp(self.mode, 'current')
                error('current injection is not yet implemented');
            end
            
            if self.stackOutput
                Inj = self.IDmeas(1:end/2) + 1i*self.IDmeas(end/2+1:end);
            else
                Inj = self.IDmeas;
            end
            if self.recalc || ~isempty(est.zeta)
                self.S = self.intS{1}/tzeta(1);
                for ii=2:length(self.intS)
                    self.S = self.S + self.intS{ii}/tzeta(ii);
                end
                if strcmp(self.mode, 'potential')
                    tempM = self.intM*diag(1./tzeta);
                    self.S = [self.S zeros(size(self.S,1), size(self.C,2));...
                              -self.C'*tempM' self.C'*self.C];
                    ninj = numel(Inj)/(self.fmesh.nel);
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(Inj, self.fmesh.nel, ninj);
                    self.b(1:self.fmesh.ng,:) = tempM*InjMat;
                    self.b(end-self.fmesh.nel+2:end,:) = ...%self.b(end-self.fmesh.nel+2:end,:)...
                                        - self.C'*(diag(self.intB./tzeta)*InjMat);
                elseif strcmp(self.mode, 'current')%current injection is not yet implemented!
                    C2 = -self.C'*diag(1./tzeta)*self.intM';
                    self.S = [self.S C2'; C2 -self.C'*diag(self.intB./tzeta)*self.C];
                    ninj = numel(Inj)/self.fmesh.nel;
                    self.b = zeros(size(self.S,1), ninj);
                    InjMat = reshape(Inj, self.fmesh.nel, ninj);
                    self.b(end-ninj+2:end,:) = self.C'*InjMat;
                else
                    error(['Unrecognized solver mode: ' self.mode]);
                end
                self.recalc = 0;
            end
            
            Cnew = [sparse(self.fmesh.ng, self.fmesh.ng) sparse(self.fmesh.ng, self.fmesh.nel-1); -self.C'*Cc' sparse(self.fmesh.nel-1, self.fmesh.nel-1)];

            self.A = A0 + self.S + Cnew;

            rhs(self.fmesh.ng+1:end,:) = - self.C'*diag(elcs);
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
        
        function vec = SolveForwardVec(self, est)
            %Solve FEM with given sigma.
            %output: vec = the currents or potentials in vector format
            %              computed using FEM
            
            %if nargin == 2 %assume "sigma" is actually "coupling"
            %    coupling = sigma;
            %    sigma = self.sigma;
            %end
            
            vec = self.SolveForward(est);
            vec = vec(:);
            vec = vec(self.vincl);
            vec = vec + self.eps;%Add a given constant (scalar or vector) to the results
            if self.stackOutput
                vec = [real(vec); imag(vec)];
            end
        end %end solveForwardVec
        
        function res = OptimizationFunction(self, est)
            %Calculate the residual of the forward problem, which is to be
            %minimized (in addition to the regularization) when solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)        
           
            elval = self.SolveForwardVec(est);
            res = 0.5*sum((self.Ln*(abs(self.Dmeas - elval))).^2);
        end
        
        
        function [Hess, grad] = GetHessAndGrad(self, est)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)
            
            est2 = DLSS_est();
            if ~isempty(est.coupl)
                if size(est.coupl,1) > 1
                    est2.coupl = est.coupl;
                else
                    est2.coupl = est.coupl*ones(size(self.fmesh.itofc,2),1);
                end
                self.coupl = est2.coupl;
            end
            if ~isempty(est.sigma)
                if size(est.sigma,1) > 1
                    est2.sigma = est.sigma;
                else
                    est2.sigma = est.sigma*ones(size(self.fmesh.itof,2),1);
                end
                self.sigma = est2.sigma;
            end
            if ~isempty(est.zeta)
                if size(est.zeta,1) > 1
                    est2.zeta = est.zeta;
                else
                    est2.zeta = est.zeta*ones(self.fmesh.nel,1);
                end
                self.zeta = est2.zeta;
            end
            fres = self.SolveForwardVec(est2);
            diff = (fres - self.Dmeas);
            Hess = DLSS_Hess();
            grad = DLSS_est();
            if ~isempty(est.coupl)
                Jc = self.FullJacobianCoupl(est2);
                Hess.coupl = Jc'*self.InvGamma_n*Jc;
                grad.coupl = Jc'*self.InvGamma_n*diff;
                if size(est.coupl,1) == 1
                    Hess.coupl = sum(sum(Hess.coupl));
                    grad.coupl = sum(grad.coupl);
                end
            end
            if ~isempty(est.sigma)
                Js = self.FullJacobianSigma(est2);
                Hess.sigma = Js'*self.InvGamma_n*Js;
                grad.sigma = Js'*self.InvGamma_n*diff;
                if size(est.sigma,1) == 1
                    Hess.sigma = sum(sum(Hess.sigma));
                    grad.sigma = sum(grad.sigma);
                end
            end
            if ~isempty(est.zeta)
                Jz = self.FullJacobianZeta(est2);
                Hess.zeta = Jz'*self.InvGamma_n*Jz;
                grad.zeta = Jz'*self.InvGamma_n*diff;
                if size(est.zeta,1) == 1
                    Hess.zeta = sum(sum(Hess.zeta));
                    grad.zeta = sum(grad.zeta);
                end
            end

            if ~isempty(est.sigma) && ~isempty(est.coupl)
                Hess.sc = Js'*self.InvGamma_n*Jc;
                if size(est.sigma) == 1
                    Hess.sc = sum(Hess.sc,1);
                end
                if size(est.coupl) == 1
                    Hess.sc = sum(Hess.sc,2);
                end
            end
            if ~isempty(est.zeta) && ~isempty(est.coupl)
                Hess.cz = Jc'*self.InvGamma_n*Jz;
                if size(est.coupl) == 1
                    Hess.cz = sum(Hess.cz,1);
                end
                if size(est.zeta) == 1
                    Hess.cz = sum(Hess.cz,2);
                end
            end
            if ~isempty(est.zeta) && ~isempty(est.sigma)
                Hess.sz = Js'*self.InvGamma_n*Jz;
                if size(est.sigma) == 1
                    Hess.sz = sum(Hess.sz,1);
                end
                if size(est.zeta) == 1
                    Hess.sz = sum(Hess.sz,2);
                end
            end
            
        end

        function [Hess, grad, Hessmix] = GetHessAndGradSigma(self, coupl, sigma)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)
            error('This should not be used anymore?');
            
            self.sigma = sigma;
            J = self.FullJacobian(sigma);
            Hess = J'*self.InvGamma_n*J;%First order approximation of the forward problem
            fres = self.SolveForwardVec(sigma, coupl);
            grad = J'*self.InvGamma_n*(fres - self.Dmeas);
            Jc = self.FullJacobianCoupl(coupl);
            Hessmix = J'*self.InvGamma_n*Jc;
            
        end

        function [Hess, grad] = GetHessAndGradCoupl(self, coupling)
            %Calculate the Hess-matrix and gradient used in solving the
            %inverse problem. This function is called by the inverse
            %problem solver (e.g. SolverLinesearch.m)
            error('This should not be used anymore?');
            
            J = self.FullJacobianCoupl(coupling);
            Hess = J'*self.InvGamma_n*J;%First order approximation of the forward problem
            fres = self.SolveForwardVec(coupling);
            grad = J'*self.InvGamma_n*(fres - self.Dmeas);
            
        end
        
        function J = FullJacobianCoupl(self, est)
            %Calculate the Jacobian of the forward problem at sigma            
            
            J = self.JacobianCoupl(est);
            if length(self.scales) == 2%The scaling done before solving the FEM has to be taken into account here.
                J(:,1:end/2) = J(:,1:end/2)*self.scales(1);
                J(:,end/2+1:end) = J(:,end/2+1:end)*self.scales(2);
            else
                J = J*diag(self.scales);%The scaling done before solving the FEM has to be taken into account here.
            end
            
        end
        
        function J = JacobianCoupl(self, est)
            %Calculate the Jacobian of the forward problem at sigma with
            %injection frequency omega. For info on mode and meshtype see
            %function SolveForward
            
            elval = self.SolveForward(est);%This elval may be inaccurate since noscf=1, but the values are not used anywhere
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
             
        function J = FullJacobianSigma(self, est)
            J = self.FullJacobian(est.sigma);
        end

        function J = FullJacobianZeta(self, est)
            %Calculate the Jacobian of the forward problem at sigma            
            
            J = self.JacobianZeta(est);
            if length(self.scales) == 2%The scaling done before solving the FEM has to be taken into account here.
                J(:,1:end/2) = J(:,1:end/2)*self.scales(1);
                J(:,end/2+1:end) = J(:,end/2+1:end)*self.scales(2);
            else
                J = J*diag(self.scales);%The scaling done before solving the FEM has to be taken into account here.
            end
            
        end

        function J = JacobianZeta(self, est)
            %Calculate the Jacobian of the forward problem at sigma with
            %injection frequency omega. For info on mode and meshtype see
            %function SolveForward
            
            elval = self.SolveForward(est);%This elval may be inaccurate since noscf=1, but the values are not used anywhere
            %if exist('Jacobian3dECTmx.mexw64')
            if 0==1%Using the mex function is currently disabled for easier debugging
                Jleft = -full(self.QC/self.A);%Using the mex function requires some pre-processing
                tempStruct.C = self.QC;
                tempStruct.Electrode = Iel;
                tempStruct.Pot = self.Pot(1:self.ng, :);
                tempStruct.PotF = self.Pot;
                Js = Jacobian3dECTmx(self.Ai, self.Av, Jleft, tempStruct);
            else
                Js = self.Jacobian3dECTZeta(elval, est);
            end
            J = Js(self.vincl,:);%some post processing of the Jacobian
            %J = self.fmesh.JacobianFtoIc(Js);
            
            
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
        
        function Js = Jacobian3dECTZeta(self, elval, est)
            % Computes the Jacobian J = d(measurements) / d(zeta)
            % Original CEM-ECT version written by K. Karhunen 03.06.2013 is implementation
            % from the earlier versions by V. Rimpiläinen(?) and L. M. Heikkinen.
            % Current version adapted for use by P. Kuusela 15.6.2023
            b = length(self.zeta);
            c = size(elval,1); 
            d = size(self.Pot,2); 
            
            Uel = reshape(self.IDmeas, c, d);
            Jleft  = self.QC/self.A;

            Jright = self.Pot;

            Js = zeros(c*d,b);

            for ii=1:b
              rzeta = zeros(size(self.zeta));
              rzeta(ii) = 1;
              tempM = self.intM*diag(rzeta);
              tempS = [self.intS{ii} zeros(size(self.intS{ii},1), size(self.C,2));...
                              -self.C'*tempM' zeros(size(self.C,2))];
              ninj = numel(Uel)/(self.fmesh.nel);
              tb = zeros(size(tempS,1), ninj);
              %InjMat = reshape(self.IDmeas, self.fmesh.nel, ninj);
              tb(1:self.fmesh.ng,:) = tempM*Uel;
              tb(end-self.fmesh.nel+2:end,:) = ...
                                        - self.C'*(diag(self.intB.*rzeta)*Uel);
              Jtemp   = -1/est.zeta(ii)^2*Jleft*(tb - tempS*Jright);
              Js(:,ii) = Jtemp(:);
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
classdef ForwardMeshDLSS < ForwardMesh2D1st
%This is a forward-mesh class for double layer sensing skin

    
    properties
        

        ng1
        itofc
        ci
        offset
        eln
        elng
        
    end
    
    methods
        
        function obj = ForwardMeshDLSS(g, H, E, ng1, ci, eln, offset)
            
            obj@ForwardMesh2D1st(g, H, E);
            obj.ng1 = ng1;
            obj.ci = ci;
            obj.itofc = 1;
            obj.eln = eln;
            if nargin < 7
                obj.offset = [1; 0];
            else
                obj.offset = offset;
            end
            
            
            
        end
        
        function [A, rhs, C2, elcs] = SigmadPhiidPhij(self, sigmaorig, couplingorig, Uel)
           %This function returns matrix A containing in its element (i,j)
           %the integral of sigma*grad(phi_i) dot grad(phi_j).
            
            A = SigmadPhiidPhij@ForwardMesh2D1st(self, sigmaorig);

            coupling = self.itofc*couplingorig;


            k = 1;  
            Arow = zeros(6*self.nH,3);
            Acol = zeros(6*self.nH,3);
            Aval = zeros(6*self.nH,3);

            k3 = 1;  
            Crow = zeros(self.nH,1);
            Ccol = zeros(self.nH,1);
            Cval = zeros(self.nH,1);

            rhsrow = zeros(self.nH,size(Uel,2));
            rhscol = zeros(self.nH,size(Uel,2));
            rhsval = zeros(self.nH,size(Uel,2));

            % Gauss quadrature points and weights
            ip=[0 0;
                1 0;
                0 1;
                0.5 0;
                0.5 0.5;
                0 0.5;
                1/3 1/3];
            a = 3/120;
            b = 8/120;
            c = 27/120;
            w = [a; a; a; b; b; b; c];

            % difference matrix of the (linear) basis functions
            L=[-1 1 0;-1 0 1];
            nInj = size(Uel,2);
            k2 = 1;
            for ii=1:self.nH
              % Go through all tetrahedra
              ind = self.H(ii,:);
              gg = self.g(ind,:);
              cc = coupling(ind);
              int = triangcphiiphij(gg,cc,ip,L,w);
              id = ind(:);
              id = [id id id];
              id2 = [ind;ind;ind];
              if self.eln(ii) == 0
                  Arow(k:k+2,:) = id;
                  Acol(k:k+2,:) = id2; 
                  Aval(k:k+2,:) = int;
                  Arow(k+3:k+5,:) = id;
                  Acol(k+3:k+5,:) = self.ci(id2);
                  Aval(k+3:k+5,:) = -int;
                  k = k + 6;
              else
                  Arow(k:k+2,:) = id;
                  Acol(k:k+2,:) = id2; 
                  Aval(k:k+2,:) = int;
                  k = k + 3;
                  rhsrow(k2:k2+2,:) = repmat(ind',1,nInj);
                  rhscol(k2:k2+2,:) = repmat(1:nInj,3,1);
                  rhsval(k2:k2+2,:) = sum(int)'*Uel(self.eln(ii),:);
                  k2 = k2 + 3;
                  Crow(k3:k3+2) = id(:,1);
                  Ccol(k3:k3+2) = self.eln(ii);
                  Cval(k3:k3+2) = sum(int)';
                  k3 = k3 + 3;
              end
            end  
            
            Aval(k:end,:) = [];
            Arow(k:end,:) = [];
            Acol(k:end,:) = [];
            rhsval(k2:end,:) = [];
            rhsrow(k2:end,:) = [];
            rhscol(k2:end,:) = [];

            
            A2 = sparse(Arow,Acol,Aval,self.ng+self.nel-1,self.ng+self.nel-1);
            A = A + A2;
            rhs = sparse(rhsrow,rhscol,rhsval,self.ng+self.nel-1,size(Uel,2));
            Cval(Ccol==0) = [];
            Crow(Crow==0) = [];
            Ccol(Ccol==0) = [];
            C2 = sparse(Crow,Ccol,Cval,self.ng,self.nel);

            %HOXHOXHO TODO this part only works with excitation pattern eye
            elcs = sum(rhs);

            
        end
        
        
        
        function [Ai, Av] = GradientMatrixCoupling(self)
            %not actually gradients, but whatevs...
            %Computes all the basis function gradients of the 1st order
            %mesh to be used in calculating the Jacobian.
            %Based on earlier UEF codes. Modified to fit EITFEM-class by
            %Petri Kuusela 9.12.2021.
            Arow = zeros(100,self.gdim+1);%guesstimate of required size
            Acol = zeros(100,self.gdim+1);
            Aval = zeros(100,self.gdim+1);

            Ai = cell(self.ng,1);
            Av = cell(self.ng,1);
            
            %constants needed in integration
            ip=[0 0;
                1 0;
                0 1;
                0.5 0;
                0.5 0.5;
                0 0.5;
                1/3 1/3];
            a = 3/120;
            b = 8/120;
            c = 27/120;
            w = [a; a; a; b; b; b; c];
            L=[-ones(self.gdim,1) eye(self.gdim)];
            
            % compute elements
            for jj=1:self.ng

              El = self.EC(jj,self.EC(jj,:)>0);

              k = 1;
              for ii=1:length(El)
                ind = self.H(El(ii),:); % Indices of the element
                gg=self.g(ind,:);
                cc = (ind == jj);
                cc = cc';
                if ~any(cc)
                    error('something wrong with element connections (EC)');
                end

                ind = repmat(ind,self.gdim+1,1);
                id = ind';
                
                int = triangcphiiphij(gg,cc,ip,L,w);
                % temporary storage
                Arow(k:k+2,:) = id;
                Acol(k:k+2,:) = ind;  
                Aval(k:k+2,:) = int;
                k = k + 3;
                if ~any(any(self.ci(ind) == 0)) && ~any(any(self.ci(ind)<8))
                    Arow(k:k+2,:) = id;
                    Acol(k:k+2,:) = self.ci(ind);
                    Aval(k:k+2,:) = -int;
                    Arow(k+3:k+5,:) = self.ci(id);
                    Acol(k+3:k+5,:) = ind;
                    Aval(k+3:k+5,:) = -int;
                    k = k + 6;
                end     
              end     
              I = 1:(k-1);
              
              tArow = Arow(I,:);
              tAcol = Acol(I,:);
              tAval = Aval(I,:);
              %tAval(tAcol == 0) = [];
              %tArow(tAcol == 0) = [];
              %tAcol(tAcol == 0) = [];
              %tAval(tArow == 0) = [];
              %tAcol(tArow == 0) = [];
              %tArow(tArow == 0) = [];
              S = sparse(tArow,tAcol,tAval,self.ng,self.ng);
              [I,~] = find(S);

              LL = unique(I);    
              Ai{jj} = LL;             % symmetric assumption => no columns needed
              Av{jj} = full(S(LL,LL));
            end
            
            


        end

        function Jr = RhsJacobianTerms(self)

            Jr = cell(self.ng,1);
            
            %constants needed in integration
            ip=[0 0;
                1 0;
                0 1;
                0.5 0;
                0.5 0.5;
                0 0.5;
                1/3 1/3];
            a = 3/120;
            b = 8/120;
            c = 27/120;
            w = [a; a; a; b; b; b; c];
            L=[-ones(self.gdim,1) eye(self.gdim)];
            
            %go through the elements
            for iH = 1:self.nH
                if self.eln(iH) == 0
                    continue;
                end
                ind = self.H(iH,:);
                gg = self.g(ind,:);
                self.elng(ind) = self.eln(iH);
                for indii = 1:length(ind)
                    ii = ind(indii);
                    cc = (ind == ii);
                    cc = cc';
                    if ~any(cc)
                        error('something went wrong!');
                    end
    
                    ind2 = repmat(ind,self.gdim+1,1);
                    id = ind2';
                    
                    int = triangcphiiphij(gg,cc,ip,L,w);
                    
                    if isempty(Jr{ii})
                        Jr{ii} = sparse(self.ng+self.nel-1,1);
                    end
                    Jr{ii}(id) = Jr{ii}(id) + sum(int)';
                end

            end
 
        end

        function Jr = RhsJacobianTermsZeta(self)

            Jr = cell(self.nel,1);
            
            %constants needed in integration
            ip=[0 0;
                1 0;
                0 1;
                0.5 0;
                0.5 0.5;
                0 0.5;
                1/3 1/3];
            a = 3/120;
            b = 8/120;
            c = 27/120;
            w = [a; a; a; b; b; b; c];
            L=[-ones(self.gdim,1) eye(self.gdim)];
            
            %go through the elements
            for iH = 1:self.nH
                if self.eln(iH) == 0
                    continue;
                end
                ind = self.H(iH,:);
                gg = self.g(ind,:);
                self.elng(ind) = self.eln(iH);
                for indii = 1:length(ind)
                    ii = ind(indii);
                    cc = (ind == ii);
                    cc = cc';
                    if ~any(cc)
                        error('something went wrong!');
                    end
    
                    ind2 = repmat(ind,self.gdim+1,1);
                    id = ind2';
                    
                    int = triangcphiiphij(gg,cc,ip,L,w);
                    
                    if isempty(Jr{ii})
                        Jr{ii} = sparse(self.ng+self.nel-1,1);
                    end
                    Jr{ii}(id) = Jr{ii}(id) + sum(int)';
                end

            end
 
        end

        function SetInverseMeshCoupling(self, imesh)
            
            
            if size(imesh.g,2) == 2%Check the dimension of the inverse mesh
                itof1 = interpolatematrix2d_new(imesh.H, imesh.g, self.g(1:self.ng1,1:2));
                itof2 = interpolatematrix2d_new(imesh.H, imesh.g+self.offset, self.g(self.ng1+1:end,1:2));
                self.itofc = [itof1;itof2];
            elseif size(imesh.g,2) == 3
                error('The support for 3D inverse mesh has not yet been added');
            else
                error(['Unfit second dimension of ginv: size(ginv,2) = ' num2str(size(ginv,2))]);
            end
            %self.nginv = size(imesh.g,1);
            
        end
        

        function s1 = ItoFc(self, coupl)
            s1 = self.itofc*coupl;
        end
        function JI = JacobianFtoIc(self, J)
            JI = J*self.itofc;
        end
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
end
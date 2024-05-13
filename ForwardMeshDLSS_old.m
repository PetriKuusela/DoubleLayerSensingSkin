classdef ForwardMeshDLSS < ForwardMesh2D1st
%This is a forward-mesh class for double layer sensing skin


    
    methods
        
        function obj = ForwardMeshDLSS(g, H, E)
            
            obj@ForwardMesh2D1st(g, H, E);
            
        end
        
        function A = SigmadPhiidPhij(self, sigma, coupling)
           %This function returns matrix A containing in its element (i,j)
           %the integral of sigma*grad(phi_i) dot grad(phi_j).
            
            sigma = self.itof*sigma;
            coupling = self.itof*coupling;

            k = 1;  
            Arow = zeros(12*self.nH,3);
            Acol = zeros(12*self.nH,3);
            Aval = zeros(12*self.nH,3);

            % Gauss quadrature points and weights
            ip = [0.5 0; 0.5 0.5; 0 0.5];
            ip2=[0 0;
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
            for ii=1:self.nH
              % Go through all tetrahedra
              ind = self.H(ii,:);
              gg = self.g(ind,:);
              ss = sigma(ind);
              cc = coupling(ind);
              int = triangLinSigma(gg,ss,ip,L);
              int2 = triangcphiiphij(gg,cc,ip2,L,w);
              id = ind(:);
              id = [id id id];
              ind = [ind;ind;ind];
              Arow(k:k+2,:) = id;
              Acol(k:k+2,:) = ind; 
              Aval(k:k+2,:) = int+int2;
              Arow(k+3:k+5,:) = id+self.ng;
              Acol(k+3:k+5,:) = ind+self.ng;
              Aval(k+3:k+5,:) = int+int2;
              Arow(k+6:k+8,:) = id+self.ng;
              Acol(k+6:k+8,:) = ind;
              Aval(k+6:k+8,:) = -int2;
              Arow(k+9:k+11,:) = id;
              Acol(k+9:k+11,:) = ind+self.ng;
              Aval(k+9:k+11,:) = -int2;

              k = k + 12;
            end  
            
            A = sparse(Arow,Acol,Aval,2*self.ng+2*self.nel-1,2*self.ng+2*self.nel-1);

            
        end
        
        
        
        function [Ai, Av] = GradientMatrix(self)
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
                Arow(k+3:k+5,:) = id+self.ng;
                Acol(k+3:k+5,:) = ind+self.ng;
                Aval(k+3:k+5,:) = int;
                Arow(k+6:k+8,:) = id+self.ng;
                Acol(k+6:k+8,:) = ind;
                Aval(k+6:k+8,:) = -int;
                Arow(k+9:k+11,:) = id;
                Acol(k+9:k+11,:) = ind+self.ng;
                Aval(k+9:k+11,:) = -int;

                k = k + 12;      
              end     
              I = 1:(k-1);

              S = sparse(Arow(I,:),Acol(I,:),Aval(I,:),2*self.ng,2*self.ng);
              [I,~] = find(S);

              LL = unique(I);    
              Ai{jj} = LL;             % symmetric assumption => no columns needed
              Av{jj} = full(S(LL,LL));
            end
        end
        
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
end
classdef ForwardMeshDLSS_varmuus < ForwardMesh2D1st
%This is a forward-mesh class for double layer sensing skin

    
    properties
        
        elH    %cell array containing triangulation for electrode areas
        elg    %cell array for node points in electrode areas
        g1
        g2
        els    %indicator for does electrode belong to surface 1 or 2
        nHtot
        nHel
        ng1
        ngtot
        itofs
        itof2
        itofs2
        EC1
        EC2
        gel1
        gel2
        
    end
    
    methods
        
        function obj = ForwardMeshDLSS_varmuus(g, H, E, elH, elg1, elg2, els)
            
            obj@ForwardMesh2D1st(g, H, E);
            obj.elH = elH;
            obj.els = els;
            obj.itof2 = 1;
            
            obj.nHtot = obj.nH;
            obj.g1 = elg1;
            obj.g2 = elg2;
            obj.ng1 = size(elg1,1);
            obj.ngtot = size(elg1,1) + size(elg2,1);
            for ie = 1:length(elH)
                obj.nHtot = obj.nHtot + size(elH{ie},1);
            end
            obj.nHel = obj.nHtot - obj.nH;
            obj.SetEC();%Do this second time now that we have set obj.ng1 and such...

            
        end
        
        function [A, rhs] = SigmadPhiidPhij(self, sigmaorig, couplingorig, Uel)
           %This function returns matrix A containing in its element (i,j)
           %the integral of sigma*grad(phi_i) dot grad(phi_j).
            
            sigma = self.itofs*sigmaorig;
            coupling = self.itof*couplingorig;
            sigma2 = self.itofs2*sigmaorig;
            coupling2 = self.itof2*couplingorig;


            k = 1;  
            Arow = zeros(12*self.nH+3*self.nHel,3);
            Acol = zeros(12*self.nH+3*self.nHel,3);
            Aval = zeros(12*self.nH+3*self.nHel,3);

            rhsrow = zeros(3*self.nHel,size(Uel,2));
            rhscol = zeros(3*self.nHel,size(Uel,2));
            rhsval = zeros(3*self.nHel,size(Uel,2));

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
              ss2 = sigma2(ind+self.ng1);
              cc = coupling(ind);
              int = triangLinSigma(gg,ss,ip,L);
              ints2 = triangLinSigma(gg,ss2,ip,L);
              int2 = triangcphiiphij(gg,cc,ip2,L,w);
              id = ind(:);
              id = [id id id];
              ind = [ind;ind;ind];
              Arow(k:k+2,:) = id;
              Acol(k:k+2,:) = ind; 
              Aval(k:k+2,:) = int+int2;
              Arow(k+3:k+5,:) = id+self.ng1;
              Acol(k+3:k+5,:) = ind+self.ng1;
              Aval(k+3:k+5,:) = ints2+int2;
              Arow(k+6:k+8,:) = id+self.ng1;
              Acol(k+6:k+8,:) = ind;
              Aval(k+6:k+8,:) = -int2;
              Arow(k+9:k+11,:) = id;
              Acol(k+9:k+11,:) = ind+self.ng1;
              Aval(k+9:k+11,:) = -int2;

              k = k + 12;
            end  
            
            %Electrodes where layer 1 has nodes
            nInj = size(Uel,2);
            k2 = 1;
            for ie = 1:self.nel
                if self.els(ie) ~= 1
                    continue
                end

                for ii = 1:size(self.elH{ie},1)
                    % Go through all tetrahedra
                  ind = self.elH{ie}(ii,:);
                  gg = self.g1(ind,:);
                  ss = sigma(ind);
                  cc = coupling(ind);
                  int = triangLinSigma(gg,ss,ip,L);
                  int2 = triangcphiiphij(gg,cc,ip2,L,w);
                  rhsrow(k2:k2+2,:) = repmat(ind',1,nInj);
                  rhscol(k2:k2+2,:) = repmat(1:nInj,3,1);
                  rhsval(k2:k2+2,:) = sum(int2)'*Uel(ie,:);
                  k2 = k2 + 3;
                  id = ind(:);
                  id = [id id id];
                  ind = [ind;ind;ind];
                  Arow(k:k+2,:) = id;
                  Acol(k:k+2,:) = ind; 
                  Aval(k:k+2,:) = int+int2;
                  k = k + 3;
                    
                end
                
            end
            
            %Electrodes where layer 2 has nodes
            for ie = 1:self.nel
                if self.els(ie) ~= 2
                    continue
                end

                for ii = 1:size(self.elH{ie},1)
                    % Go through all tetrahedra
                  ind = self.elH{ie}(ii,:);
                  gg = self.g2(ind,:);
                  ss = sigma2(ind+self.ng1);
                  cc = coupling2(ind);
                  int = triangLinSigma(gg,ss,ip,L);
                  int2 = triangcphiiphij(gg,cc,ip2,L,w);
                  rhsrow(k2:k2+2,:) = repmat(ind',1,nInj)+self.ng1;
                  rhscol(k2:k2+2,:) = repmat(1:nInj,3,1);
                  rhsval(k2:k2+2,:) = sum(int2)'*Uel(ie,:);
                  k2 = k2 + 3;
                  id = ind(:);
                  id = [id id id];
                  ind = [ind;ind;ind];
                  Arow(k:k+2,:) = id+self.ng1;
                  Acol(k:k+2,:) = ind+self.ng1; 
                  Aval(k:k+2,:) = int+int2;
                  k = k + 3;
                    
                end
                
            end
            
            A = sparse(Arow,Acol,Aval,self.ngtot+self.nel-1,self.ngtot+self.nel-1);
            rhs = sparse(rhsrow,rhscol,rhsval,self.ngtot+self.nel-1,size(Uel,2));
            
        end


        function [S, M, B] = EITElectrodeTerms(self)
            %This function returns matrices S and M, and vector B.
            %S is a cell array containing integrals phi_i*phi_j along
            %each electrode surfaces (one electrode per cell)
            %M contains integral phi_i on surface of electrode j
            %B contains the electrode areas in its elements.

            B = zeros(self.nel,1);

            M = zeros(self.ngtot, self.nel);

            S = cell(self.nel,1);
            % Loop through electrodes
            for ii=1:self.nel
              if self.els(ii) == 1%here the selection is in a sense "reversed" from when we consider to which layer the elements belong to
                  indadd = self.ng1;
              else
                  indadd = 0;
              end
              spos = 1;
              faces = self.E{ii};
              len = 2*size(self.E{ii},1);
              intS = zeros(len,2);
              rowS = zeros(len,2);
              colS = zeros(len,2);
              % Loop through faces on electrode ii
              for jj = 1:size(faces,1)
                ind = faces(jj,:); % face nodes

                gg = self.g(ind,:);
                rcidx = [ind;ind];

                bb1 = phii1D(gg);
                bb2 = phiiphij1D(gg);

                intS(spos:spos+1,:) = bb2;
                rowS(spos:spos+1,:) = rcidx.' + indadd; 
                colS(spos:spos+1,:) = rcidx + indadd;   

                M(ind+indadd,ii) = M(ind+indadd,ii) + bb1;
                %intM(mpos,:) = bb1.';
                %colM(mpos,:) = ind; 

                B(ii,:)  = B(ii,:) + ElectrodeArea1D(gg);%check is this sign ok?

                spos = spos + 2;
              end
              S{ii} = sparse(rowS,colS,intS,self.ngtot,self.ngtot);
            end
            
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
                Arow(k+3:k+5,:) = id+self.ng1;
                Acol(k+3:k+5,:) = ind+self.ng1;
                Aval(k+3:k+5,:) = int;
                Arow(k+6:k+8,:) = id+self.ng1;
                Acol(k+6:k+8,:) = ind;
                Aval(k+6:k+8,:) = -int;
                Arow(k+9:k+11,:) = id;
                Acol(k+9:k+11,:) = ind+self.ng1;
                Aval(k+9:k+11,:) = -int;

                k = k + 12;      
              end     
              I = 1:(k-1);

              S = sparse(Arow(I,:),Acol(I,:),Aval(I,:),self.ngtot,self.ngtot);
              [I,~] = find(S);

              LL = unique(I);    
              Ai{jj} = LL;             % symmetric assumption => no columns needed
              Av{jj} = full(S(LL,LL));
            end
            
            


        end

        function [b1, b2] = JacobianTerms(self)

            b1 = cell(self.ng1,2);
            b2 = cell(self.ngtot-self.ng1,2);
            
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
            
            %go through the electrodes
            for ie = 1:self.nel

                for jj = 1:size(self.elH{ie},1)
                    ind = self.elH{ie}(jj,:); % Indices of the element
                    if self.els(ie) == 1
                        gg=self.g1(ind,:);
                    else
                        gg=self.g2(ind,:);
                    end
                    for indii = 1:length(ind)
                        ii = ind(indii);
                        cc = (ind == ii);
                        cc = cc';
                        if ~any(cc)
                            error('something wrong with element connections (EC)');
                        end
        
                        %ind2 = repmat(ind,self.gdim+1,1);
                        %id2 = ind2';
                        id = ind';
                        
                        int = triangcphiiphij(gg,cc,ip,L,w);
                        
                        if self.els(ie) == 1
                            if isempty(b1{ii,1})
                                b1{ii,1} = sparse(self.ngtot+self.nel-1,1);
                                b1{ii,2} = sparse(self.ngtot+self.nel-1,self.ngtot+self.nel-1);
                            end
                            b1{ii,1}(id) = b1{ii,1}(id) + sum(int)';
                            b1{ii,2}(id, id) = b1{ii,2}(id, id) + int;
                        elseif self.els(ie) == 2
                            if isempty(b2{ii,1})
                                b2{ii,1} = sparse(self.ngtot+self.nel-1,1);
                                b2{ii,2} = sparse(self.ngtot+self.nel-1,self.ngtot+self.nel-1);
                            end
                            b2{ii,1}(id+self.ng1) = b2{ii,1}(id+self.ng1) + sum(int)';
                            b2{ii,2}(id, id) = b2{ii,2}(id, id) + int;
                        end
                    end
                end

            end


            
        end

        function [Ai, Av] = GradientMatrixForSigma(self)
            %Computes all the basis function gradients of the 1st order
            %mesh to be used in calculating the Jacobian.
            %Based on earlier UEF codes. Modified to fit EITFEM-class by
            %Petri Kuusela 9.12.2021.
            ar = zeros(self.ng*12,self.gdim+1);
            ac = zeros(self.ng*12,self.gdim+1);
            av = zeros(self.ng*12,self.gdim+1);

            Ai = cell(self.ngtot,1);
            Av = cell(self.ngtot,1);

            % compute gradients
            for jj=1:self.ng

              El = self.EC(jj,self.EC(jj,:)>0);

              rid = 1;
              for ii=1:length(El)
                ind = self.H(El(ii),:); % Indices of the element
                gg=self.g(ind,:);

                idc = repmat(ind,self.gdim+1,1);
                idr = idc';
                
                L=[-ones(self.gdim,1) eye(self.gdim)];
                Jt=L*gg;
                dJt=abs(det(Jt)); % Tetrahedra volume
                G=Jt\L; % Gradients of each basis function
                GdJt=G'*G*dJt;

                if self.gdim == 3
                    int=1/24*GdJt;
                elseif self.gdim == 2
                    int=1/6*GdJt;
                end

                % temporary storage
                ar(rid:rid+self.gdim,:) = idr;
                ac(rid:rid+self.gdim,:) = idc;
                av(rid:rid+self.gdim,:) = int;
                ar(rid+self.gdim+1:rid+2*self.gdim+1,:) = idr + self.ng1;
                ac(rid+self.gdim+1:rid+2*self.gdim+1,:) = idc + self.ng1;
                av(rid+self.gdim+1:rid+2*self.gdim+1,:) = int;
                rid = rid + 2 + 2*self.gdim;      
              end  
              ridnoelectrodes = rid;
              if self.gel1(jj) > 0
                  El = self.EC1(jj,self.EC1(jj,:)>0);
                  for ii=1:length(El)
                    ind = self.elH{self.gel1(jj)}(El(ii),:); % Indices of the element
                    gg=self.g1(ind,:);
    
                    idc = repmat(ind,self.gdim+1,1);
                    idr = idc';
                    
                    L=[-ones(self.gdim,1) eye(self.gdim)];
                    Jt=L*gg;
                    dJt=abs(det(Jt)); % Tetrahedra volume
                    G=Jt\L; % Gradients of each basis function
                    GdJt=G'*G*dJt;
    
                    if self.gdim == 3
                        int=1/24*GdJt;
                    elseif self.gdim == 2
                        int=1/6*GdJt;
                    end
    
                    % temporary storage
                    ar(rid:rid+self.gdim,:) = idr;
                    ac(rid:rid+self.gdim,:) = idc;
                    av(rid:rid+self.gdim,:) = int;
                    rid = rid + 1 + self.gdim;      
                  end  
              end
              if self.gel2(jj) > 0
                  El = self.EC2(jj,self.EC2(jj,:)>0);
                  for ii=1:length(El)
                    ind = self.elH{self.gel2(jj)}(El(ii),:); % Indices of the element
                    gg=self.g2(ind,:);
    
                    idc = repmat(ind,self.gdim+1,1);
                    idr = idc';
                    
                    L=[-ones(self.gdim,1) eye(self.gdim)];
                    Jt=L*gg;
                    dJt=abs(det(Jt)); % Tetrahedra volume
                    G=Jt\L; % Gradients of each basis function
                    GdJt=G'*G*dJt;
    
                    if self.gdim == 3
                        int=1/24*GdJt;
                    elseif self.gdim == 2
                        int=1/6*GdJt;
                    end
    
                    % temporary storage
                    ar(rid:rid+self.gdim,:) = idr+self.ng1;
                    ac(rid:rid+self.gdim,:) = idc+self.ng1;
                    av(rid:rid+self.gdim,:) = int;
                    rid = rid + 1 + self.gdim;      
                  end  
              end

              I = 1:(rid-1);
              S = sparse(ar(I,:),ac(I,:),av(I,:),self.ngtot,self.ngtot);
              [I,~] = find(S);
              L = unique(I);    
              if self.gel1(jj) > 0
                Ai{jj} = L;             % symmetric assumption => no columns needed
                Av{jj} = full(S(L,L));
                I = 1:(ridnoelectrodes-1);
                S = sparse(ar(I,:),ac(I,:),av(I,:),self.ngtot,self.ngtot);
                [I,~] = find(S);
                L = unique(I);  
                Ai{jj+self.ng1} = L;             % symmetric assumption => no columns needed
                Av{jj+self.ng1} = full(S(L,L));
              elseif self.gel2(jj) > 0
                Ai{jj+self.ng1} = L;             % symmetric assumption => no columns needed
                Av{jj+self.ng1} = full(S(L,L));
                I = 1:(ridnoelectrodes-1);
                S = sparse(ar(I,:),ac(I,:),av(I,:),self.ngtot,self.ngtot);
                [I,~] = find(S);
                L = unique(I);  
                Ai{jj} = L;             % symmetric assumption => no columns needed
                Av{jj} = full(S(L,L));
              else
                Ai{jj} = L;             % symmetric assumption => no columns needed
                Av{jj} = full(S(L,L));
                Ai{jj+self.ng1} = L;             % symmetric assumption => no columns needed
                Av{jj+self.ng1} = full(S(L,L));
              end

            end
            
            %Layer1 electrode areas
            for jj=self.ng+1:self.ng1

              El = self.EC1(jj,self.EC1(jj,:)>0);

              rid = 1;
              for ii=1:length(El)
                ind = self.elH{self.gel1(jj)}(El(ii),:); % Indices of the element
                gg=self.g1(ind,:);

                idc = repmat(ind,self.gdim+1,1);
                idr = idc';
                
                L=[-ones(self.gdim,1) eye(self.gdim)];
                Jt=L*gg;
                dJt=abs(det(Jt)); % Tetrahedra volume
                G=Jt\L; % Gradients of each basis function
                GdJt=G'*G*dJt;

                if self.gdim == 3
                    int=1/24*GdJt;
                elseif self.gdim == 2
                    int=1/6*GdJt;
                end

                % temporary storage
                ar(rid:rid+self.gdim,:) = idr;
                ac(rid:rid+self.gdim,:) = idc;
                av(rid:rid+self.gdim,:) = int;
                rid = rid + 1 + self.gdim;      
              end  
              
              I = 1:(rid-1);
              S = sparse(ar(I,:),ac(I,:),av(I,:),self.ngtot,self.ngtot);
              [I,~] = find(S);
              L = unique(I);    
              Ai{jj} = L;             % symmetric assumption => no columns needed
              Av{jj} = full(S(L,L));
              
            end

            %Layer2 electrode areas
            for jj=self.ng+1:self.ngtot-self.ng1

              El = self.EC2(jj,self.EC2(jj,:)>0);

              rid = 1;
              for ii=1:length(El)
                ind = self.elH{self.gel2(jj)}(El(ii),:); % Indices of the element
                gg=self.g2(ind,:);

                idc = repmat(ind,self.gdim+1,1);
                idr = idc';
                
                L=[-ones(self.gdim,1) eye(self.gdim)];
                Jt=L*gg;
                dJt=abs(det(Jt)); % Tetrahedra volume
                G=Jt\L; % Gradients of each basis function
                GdJt=G'*G*dJt;

                if self.gdim == 3
                    int=1/24*GdJt;
                elseif self.gdim == 2
                    int=1/6*GdJt;
                end

                % temporary storage
                ar(rid:rid+self.gdim,:) = idr+self.ng1;
                ac(rid:rid+self.gdim,:) = idc+self.ng1;
                av(rid:rid+self.gdim,:) = int;
                rid = rid + 1 + self.gdim;      
              end  
              
              I = 1:(rid-1);
              S = sparse(ar(I,:),ac(I,:),av(I,:),self.ngtot,self.ngtot);
              [I,~] = find(S);
              L = unique(I);    
              Ai{jj+self.ng1} = L;             % symmetric assumption => no columns needed
              Av{jj+self.ng1} = full(S(L,L));
              
            end

        end
        


        function SetInverseMesh(self, imesh)
            
            
            %self.nginv = size(ginv,1);
            igntot = size(imesh.g1,1) + size(imesh.g2,1);
            if size(imesh.g,2) == 2%Check the dimension of the inverse mesh
                self.itof = interpolatematrix2d_new(imesh.H, imesh.g, self.g1(:,1:2));
                self.itofs = interpolatematrix2d_new(imesh.H1, imesh.g1, self.g1(:,1:2));
                self.itofs(self.ngtot,igntot) = 0;
                self.itof2 = interpolatematrix2d_new(imesh.H, imesh.g, self.g2(:,1:2));
                self.itofs2 = interpolatematrix2d_new(imesh.H2, imesh.g2, self.g2(:,1:2));
                self.itofs2 = [sparse(self.ngtot, size(imesh.g1,1)) [sparse(self.ng1, size(imesh.g2,1)); self.itofs2]];

            elseif size(imesh.g,2) == 3
                error('The support for 3D inverse mesh has not yet been added');
            else
                error(['Unfit second dimension of ginv: size(ginv,2) = ' num2str(size(ginv,2))]);
            end
            self.nginv = size(imesh.g,1);
            
        end
        

        function SetEC(self)
            %Sets self.EC, i.e. an array, for which each row i contains all
            %such element indices (i.e. rows of H) that the element contains the node i (i.e. row i of g).
            %currently the only use for self.EC is by self.GradientMatrix()
            SetEC@ForwardMesh2D1st(self);
            self.EC1 = zeros(self.ng1, 3);
            self.EC2 = zeros(self.ngtot-self.ng1, 3);
            counter1 = zeros(self.ng1, 1);
            counter2 = zeros(self.ngtot-self.ng1, 1);
            self.gel1 = zeros(self.ng1, 1);
            self.gel2 = zeros(self.ngtot-self.ng1, 1);
            for ie = 1:length(self.elH)
                for iH = 1:size(self.elH{ie},1)
                    if self.els(ie) == 1
                        for in = 1:size(self.elH{ie},2)
                            id = self.elH{ie}(iH, in);
                            self.gel1(id) = ie;
                            counter1(id) = counter1(id) + 1;
                            self.EC1(id, counter1(id)) = iH;
                        end
                    elseif self.els(ie) == 2
                        for in = 1:size(self.elH{ie},2)
                            id = self.elH{ie}(iH, in);
                            self.gel2(id) = ie;
                            counter2(id) = counter2(id) + 1;
                            self.EC2(id, counter2(id)) = iH;
                        end
                    end
                end
            end
            
        end
        
        function s1 = ItoF(self, sigma)
            s1 = [self.itof*sigma; self.itof2*sigma];
        end
        function s1 = ItoFForSigma(self, sigma)
            s1 = [self.itofs*sigma; self.itofs2*sigma];
        end
        function JI = JacobianFtoI(self, J)
            JI = J*[self.itof; self.itof2(self.ng+1:end,:)];
        end
        function JI = JacobianFtoIForSigma(self, J)
            %JI = J*[self.itofs; self.itofs2];
            JI = J*(self.itofs +  self.itofs2);
        end
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
    
end
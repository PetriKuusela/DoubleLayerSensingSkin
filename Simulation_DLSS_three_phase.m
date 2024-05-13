function [Umeas, Imeas, Umeas_i, Imeas_i, Umeas_0, Imeas_0, e] = Simulation_DLSS_three_phase(fm, sigma, coupling, z, vincl, imesh, err, mode, absmode, sigma_h, nsginv, sigma2)
%A simple simulation that will return synthetic data. meas_i refer to
%initial homogeneous solutions, whereas Umeas and Imeas are the situation
%we are actually interested in. if sigma_h = 0, reference measurements are
%not calculated.
%
%Author: Petri Kuusela 16.11.2021

%Use default values for arguments that are missing:
if nargin < 7 || isempty(err)
    err = [1e-6 1e-3 1e-4 1e-2];
end
if nargin < 8 || isempty(mode)
    mode = 'potential';
end
if nargin < 9 || isempty(absmode)
    absmode = 0;
end
if nargin < 10 || isempty(sigma_h)
    sigma_h = 0;
end
if nargin < 12 || isempty(sigma2)
    sigma2 = [];
end

    %These are the noise parameters
    noiserel = err(1);
    noiseabs = err(2);
    esystematicrel = err(3);
    esystematicabs = err(4);

    solver = EITFEM_DLSS(fm);
    solver.zeta = z;
    solver.vincl = vincl;
    solver.mode = mode;
    solver.sigmamin = 1e-9;
    if strcmp(mode, 'potential')
        solver.IDmeas = eye(length(fm.E));%Injection pattern
        solver.IDmeas= solver.IDmeas(:);
        if length(sigma) == 2*length(imesh.g)
            solver.IDmeas = [solver.IDmeas; 0*solver.IDmeas];
        end
    elseif strcmp(mode, 'current')
        skipone = 0;
        if skipone
            Imeas = [1; 0; -1; zeros(length(fm.E)-1,1)];
            Imeas = repmat(Imeas, length(fm.E),1);
            Imeas(end-length(fm.E)+1:end) = [];
            Imeas(end-length(fm.E)+1) = -1;
        else
            Imeas = eye(length(fm.E));%Create the injection pattern
            Imeas(2:end,1:end-1) = Imeas(2:end,1:end-1) - eye(length(fm.E)-1);
            Imeas(1,end) = -1;
        end
        solver.IDmeas = Imeas(:);
        if absmode == 0 && length(sigma) == 2*length(imesh.g)
            solver.IDmeas = [solver.IDmeas; 0*solver.IDmeas];
        end
        solver.IDmeas = solver.IDmeas*1e-3;%Injected currents are usually order of mA
    end
    if size(imesh.g,1) ~= size(fm.g,1)
        fm.SetInverseMeshCoupling(imesh);
    end
    if length(sigma) == 2*length(imesh.g)
        solver.stackOutput = 1;
    end

    coupl0 = zeros(length(coupling),1);
    
    if strcmp(mode, 'potential')
        Imeas_0 = solver.SolveForwardVec(sigma, coupl0);%these are the results
        Umeas_0 = solver.IDmeas;
        %add noise and error:
        esys = randn(length(Imeas_0),1)*esystematicabs*(max(Imeas_0)-min(Imeas_0));
        esysrel = randn(length(Imeas_0),1)*esystematicrel;
        Imeas_0 = Imeas_0.*(1+esysrel) + esys;
        Imeas_0 = Imeas_0.*(1+randn(length(Imeas_0),1)*noiserel) + randn(length(Imeas_0),1)*noiseabs*(max(Imeas_0(Imeas_0>0)));
        Imeas_0(Imeas_0<0 & Imeas_0>0.1*min(Imeas_0)) = 1e-6;
    elseif strcmp(mode, 'current')
        Umeas_i = solver.SolveForwardVec(sigma, coupl0);%these are the results
        Imeas_i = solver.IDmeas;
        %add noise and error:
        esys = randn(length(Umeas_0),1)*esystematicabs*(max(Umeas_0)-min(Umeas_0));
        esysrel = randn(length(Umeas_0),1)*esystematicrel;
        Umeas_0 = Umeas_0.*(1+esysrel) + esys;
        Umeas_0 = Umeas_0.*(1+randn(length(Umeas_0),1)*noiserel) + randn(length(Umeas_0),1)*noiseabs*(max(Umeas_0)-min(Umeas_0));
    end

        
    if sigma_h(1) ~= 0
        %homogeneous measurements:
        coupl_hfull = sigma_h(1)*ones(length(coupling),1);
        if strcmp(mode, 'potential')
            Imeas_i = solver.SolveForwardVec(sigma, coupl_hfull);%these are the results
            Umeas_i = solver.IDmeas;
            %add noise and error:
            if ~exist('esys')%there have been no homogeneous measurements where these have already been calculated
                esys = randn(length(Imeas_i),1)*esystematicabs*(max(Imeas_i)-min(Imeas_i));
                esysrel = randn(length(Imeas_i),1)*esystematicrel;
            end
            Imeas_i = Imeas_i.*(1+esysrel) + esys;
            Imeas_i = Imeas_i.*(1+randn(length(Imeas_i),1)*noiserel) + randn(length(Imeas_i),1)*noiseabs*(max(Imeas_i(Imeas_i>0)));
            Imeas_i(Imeas_i<0 & Imeas_i>0.1*min(Imeas_i)) = 1e-6;
        elseif strcmp(mode, 'current')
            Umeas_i = solver.SolveForwardVec(sigma, coupl_hfull);%these are the results
            Imeas_i = solver.IDmeas;
            %add noise and error:
            if ~exist('esys')%there have been no homogeneous measurements where these have already been calculated
                esys = randn(length(Umeas_i),1)*esystematicabs*(max(Umeas_i)-min(Umeas_i));
                esysrel = randn(length(Umeas_i),1)*esystematicrel;
            end
            Umeas_i = Umeas_i.*(1+esysrel) + esys;
            Umeas_i = Umeas_i.*(1+randn(length(Umeas_i),1)*noiserel) + randn(length(Umeas_i),1)*noiseabs*(max(Umeas_i)-min(Umeas_i));
        end

        
    end %end homogeneous measurements
    
    if ~isempty(sigma2)
        sigma = sigma2;
    end
        
    %The actual measurements:
    if strcmp(mode, 'potential')
        Imeas = solver.SolveForwardVec(sigma, coupling);%these are the results
        Umeas = solver.IDmeas;
        %add noise and error:
        if ~exist('esys')%there have been no homogeneous measurements where these have already been calculated
            esys = randn(length(Imeas),1)*esystematicabs*(max(Imeas)-min(Imeas));
            esysrel = randn(length(Imeas),1)*esystematicrel;
        end%If previous esys and esysrel exist use them
        Imeas = Imeas.*(1+esysrel) + esys;
        Imeas = Imeas.*(1+randn(length(Imeas),1)*noiserel) + randn(length(Imeas),1)*noiseabs*(max(Imeas(Imeas>0)));
        Imeas(Imeas<0 & Imeas>0.1*min(Imeas)) = 1e-6;
    elseif strcmp(mode, 'current')
        Umeas = solver.SolveForwardVec(sigma, coupling);%these are the results
        Imeas = solver.IDmeas;
        %add noise and error:
        if ~exist('esys')%there have been no homogeneous measurements where these have already been calculated
            esys = randn(length(Umeas),1)*esystematicabs*(max(Umeas)-min(Umeas));
            esysrel = randn(length(Umeas),1)*esystematicrel;
        end%If previous esys and esysrel exist use them
        Umeas = Umeas.*(1+esysrel) + esys;
        Umeas = Umeas.*(1+randn(length(Umeas),1)*noiserel) + randn(length(Umeas),1)*noiseabs*(max(Umeas)-min(Umeas));
    end
    
    %return also the systematic error:
    e = [esys esysrel];
    if ~exist('Umeas_i')%reference measurements were not done
        Umeas_i = 0;
        Imeas_i = 0;
    end
    
    
end
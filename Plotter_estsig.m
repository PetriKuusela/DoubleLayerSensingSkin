classdef Plotter_estsig < handle
%A Plotter class that can be assigned to inverse problem solver (e.g.
%SolverLinesearch.m) to plot the estimate in between the iterations.
    
    properties
        gs       %The array of node co-ordinates
        Hs       %The elements of the mesh
        gc
        Hc
        g1      %The nodes of the forward mesh (used to plot the electrodes)
        elfaces %A cell array containing electrode nodes referring to g1
        cm      %colormap
        scales  %the scales to multiply the estimates to get the plotted values
        fig     %A figure handle for conductivity
        fig2    %A figure handle for coupling
        fig3    %A figure handle for zetas
        plotel  %A flag: do we want to plot the electrodes?
    end
    
    methods
        
        function obj = Plotter_estsig(gs, Hs, gc, Hc, g1, elfaces)
            %Class constructor.
            %Input: g and H define the mesh we want to plot on
            %       g1 and elfaces(optional) define where the electrodes
            %       are plotted
            obj.gs = gs;
            obj.Hs = Hs;
            obj.gc = gc;
            obj.Hc = Hc;
            obj.cm = jet;
            obj.fig = figure();
            obj.fig2 = figure();
            obj.fig3 = figure();
            obj.plotel = 0;
            if nargin > 5
                obj.elfaces = elfaces;
                obj.g1 = g1;
                obj.plotel = 1;
            end
            obj.scales = 1;
        end
        
        function plot(self, est, minmaxvals)
            %The basic plot function to plot sigma. The colorbar scale can
            %be given in minmaxvals = [minval maxval].
            
            
            %First check if we try to plot the electrodes without having
            %given the necessary variables:
            if self.plotel && (isempty(self.elfaces) || isempty(self.g1))
                warning('Cannot plot electrodes since elfaces or g1 is not set');
                self.plotel = 0;
            end
            
            %Set the figure we want to plot on:
            set(0, 'CurrentFigure', self.fig);
            clf;
                
            %Apply the scales:
            %est.sigma = est.sigma.*self.scales.sigma;
            %est.coupl = est.coupl.*self.scales.coupl;

            if ~isempty(est.sigma)
                if nargin > 2
                    minval = minmaxvals(1);
                    maxval = minmaxvals(2);
                else
                    minval = min(est.sigma)-1e-9;
                    maxval = max(est.sigma)+1e-9;
                end
                
                
                fh = trisurf(self.Hs, self.gs(:,1), self.gs(:,2), est.sigma);%This is the line that does the main plotting
                view(2);
                set(fh, 'edgecolor', 'none');%no edges
                set(fh, 'FaceColor', 'interp');%smooth colors between elements
                colormap(get(fh,'parent'), self.cm);%set colormap
                if maxval > minval
                    set(get(fh,'parent'),'CLim',[minval maxval]);
                else
                    set(get(fh,'parent'),'CLim',[minval-1e-6 maxval+1e-6]);
                end
                axis equal
                colorbar;
                if self.plotel
                    hold on;
                    for iel = 1:length(self.elfaces)
                        for it = 1:length(self.elfaces{iel})
                            plot(self.g1(self.elfaces{iel}(it,:), 1), self.g1(self.elfaces{iel}(it,:), 2), 'k-', 'LineWidth', 5);
                        end
                    end
                    hold off;
                end
            
            end

            %BEGIN COUPLING PART:

            if ~isempty(est.coupl) && size(est.coupl,1) > 1
                if nargin > 2
                    minval = minmaxvals(3);
                    maxval = minmaxvals(4);
                else
                    minval = min(est.coupl)-1e-9;
                    maxval = max(est.coupl)+1e-9;
                end
    
                set(0, 'CurrentFigure', self.fig2);
                clf;
                
                fh = trisurf(self.Hc, self.gc(:,1), self.gc(:,2), est.coupl);%This is the line that does the main plotting
                view(2);
                set(fh, 'edgecolor', 'none');%no edges
                set(fh, 'FaceColor', 'interp');%smooth colors between elements
                colormap(get(fh,'parent'), self.cm);%set colormap
                if maxval > minval
                    set(get(fh,'parent'),'CLim',[minval maxval]);
                else
                    set(get(fh,'parent'),'CLim',[minval-1e-6 maxval+1e-6]);
                end
                axis('square')
                colorbar;
                if self.plotel
                    hold on;
                    for iel = 1:length(self.elfaces)
                        for it = 1:length(self.elfaces{iel})
                            plot(self.g1(self.elfaces{iel}(it,:), 1), self.g1(self.elfaces{iel}(it,:), 2), 'k-', 'LineWidth', 5);
                        end
                    end
                    hold off;
                end

            end   

            %BEGIN CONTACT IMPEDANCES:
    
            if ~isempty(est.zeta)
                
                set(0, 'CurrentFigure', self.fig3);
                clf;

                plot(est.zeta);
                

            end
                
        end
        
        
    end    
    
    
end
classdef DLSS_Hess% < handle


    properties
        sigma
        coupl
        zeta
        sc
        sz
        cz
    end


    methods

        function obj = DLSS_Hess(sigma, coupl, zeta, sc, sz, cz)
            arguments
                sigma = []
                coupl = []
                zeta = []
                sc = []
                sz = []
                cz = []
            end
            if nargin == 1 && length(sigma) == 6%in this case sigma is not sigma, but cell array containing all vals
                coupl = sigma{2};
                zeta = sigma{3};
                sc = sigma{4};
                sz = sigma{5};
                cz = sigma{6};
                sigma = sigma{1};
            end
            obj.sigma = sigma;
            obj.coupl = coupl;
            obj.zeta = zeta;
            if isempty(sc) && ~isempty(sigma) && ~isempty(coupl)
                obj.sc = zeros(size(sigma,1), size(coupl,1));
            else
                obj.sc = sc;
            end
            if isempty(sz) && ~isempty(sigma) && ~isempty(zeta)
                obj.sz = zeros(size(sigma,1), size(zeta,1));
            else
                obj.sz = sz;
            end
            if isempty(cz) && ~isempty(coupl) && ~isempty(zeta)
                obj.cz = zeros(size(coupl,1), size(zeta,1));
            else
                obj.cz = cz;
            end
        end

        function r = plus(o1, o2)
            if isa(o1, 'DLSS_Hess') && isa(o2, 'DLSS_Hess')
%                 if ~isempty(o1.sigma) && ~isempty(o2.sigma)
%                     s = o1.sigma + o2.sigma;
%                 end
%                 if ~isempty(o1.coupl) && ~isempty(o2.coupl)
%                     c = o1.coupl + o2.coupl;
%                 end
%                 if ~isempty(o1.zeta) && ~isempty(o2.zeta)
%                     c = o1.zeta + o2.zeta;
%                 end
%                 if ~isempty(o1.sc) && ~isempty(o2.sc)
%                     sc = o1.sc + o2.sc;
%                 end
%                 if ~isempty(o1.sz) && ~isempty(o2.sz)
%                     sz = o1.sz + o2.sz;
%                 end
%                 if ~isempty(o1.cz) && ~isempty(o2.cz)
%                     cz = o1.cz + o2.cz;
%                 end
                props = properties(o1);
                vals = cell(length(props),1);
                for ip = 1:length(props)
                    if ~isempty(o1.(props{ip})) && ~isempty(o2.(props{ip}))
                        vals{ip} = o1.(props{ip}) + o2.(props{ip});
                    end
                end
                r = DLSS_Hess(vals);
            elseif isa(o1, 'DLSS_Hess')
%                 if ~isempty(o1.sigma)
%                     s = o1.sigma + o2;
%                 end
%                 if ~isempty(o1.coupl)
%                     c = o1.coupl + o2;
%                 end
%                 if ~isempty(o1.zeta)
%                     z = o1.zeta + o2;
%                 end
%                 if ~isempty(o1.sc)
%                     sc = o1.sc + o2;
%                 end
%                 if ~isempty(o1.sz)
%                     sz = o1.sz + o2;
%                 end
%                 if ~isempty(o1.cz)
%                     cz = o1.cz + o2;
%                 end
%                 r = DLSS_Hess(s, c, z, sc, sz, cz);
                props = properties(o1);
                vals = cell(length(props),1);
                for ip = 1:length(props)
                    if ~isempty(o1.(props{ip}))
                        vals{ip} = o1.(props{ip}) + o2;
                    end
                end
                r = DLSS_Hess(vals);
            else
                props = properties(o2);
                vals = cell(length(props),1);
                for ip = 1:length(props)
                    if ~isempty(o2.(props{ip}))
                        vals{ip} = o2.(props{ip}) + o1;
                    end
                end
                r = DLSS_Hess(vals);
            end
        end

        function r = minus(o1, o2)
            if isa(o1, 'DLSS_Hess') && isa(o2, 'DLSS_Hess')
%                 if ~isempty(o1.sigma) && ~isempty(o2.sigma)
%                     s = o1.sigma + o2.sigma;
%                 end
%                 if ~isempty(o1.coupl) && ~isempty(o2.coupl)
%                     c = o1.coupl + o2.coupl;
%                 end
%                 if ~isempty(o1.zeta) && ~isempty(o2.zeta)
%                     c = o1.zeta + o2.zeta;
%                 end
%                 if ~isempty(o1.sc) && ~isempty(o2.sc)
%                     sc = o1.sc + o2.sc;
%                 end
%                 if ~isempty(o1.sz) && ~isempty(o2.sz)
%                     sz = o1.sz + o2.sz;
%                 end
%                 if ~isempty(o1.cz) && ~isempty(o2.cz)
%                     cz = o1.cz + o2.cz;
%                 end
                props = properties(o1);
                vals = cell(length(props),1);
                for ip = 1:length(props)
                    if ~isempty(o1.(props{ip})) && ~isempty(o2.(props{ip}))
                        vals{ip} = o1.(props{ip}) - o2.(props{ip});
                    end
                end
                r = DLSS_Hess(vals);
            elseif isa(o1, 'DLSS_Hess')
%                 if ~isempty(o1.sigma)
%                     s = o1.sigma + o2;
%                 end
%                 if ~isempty(o1.coupl)
%                     c = o1.coupl + o2;
%                 end
%                 if ~isempty(o1.zeta)
%                     z = o1.zeta + o2;
%                 end
%                 if ~isempty(o1.sc)
%                     sc = o1.sc + o2;
%                 end
%                 if ~isempty(o1.sz)
%                     sz = o1.sz + o2;
%                 end
%                 if ~isempty(o1.cz)
%                     cz = o1.cz + o2;
%                 end
%                 r = DLSS_Hess(s, c, z, sc, sz, cz);
                props = properties(o1);
                vals = cell(length(props),1);
                for ip = 1:length(props)
                    if ~isempty(o1.(props{ip}))
                        vals{ip} = o1.(props{ip}) - o2;
                    end
                end
                r = DLSS_Hess(vals);
            else
                props = properties(o2);
                vals = cell(length(props),1);
                for ip = 1:length(props)
                    if ~isempty(o2.(props{ip}))
                        vals{ip} = o1 - o2.(props{ip});
                    end
                end
                r = DLSS_Hess(vals);
            end
        end

        function r = uminus(o1)
            props = properties(o1);
            vals = cell(length(props),1);
            for ip = 1:length(props)
                if ~isempty(o1.(props{ip}))
                    vals{ip} = -o1.(props{ip});
                end
            end
            r = DLSS_Hess(vals);
        end

        function r = mtimes(o1, o2)
            if ~isa(o1, 'DLSS_Hess')
                error('Type error!');
            end
            if ~isa(o2, 'DLSS_est')
                error('Type error!')
            end
            Hess = [o1.coupl o1.sc' o1.cz; o1.sc o1.sigma o1.sz; o1.cz' o1.sz' o1.zeta];
            grad = [o2.coupl; o2.sigma; o2.zeta];
            tr = Hess*grad;
            r = DLSS_est(tr(length(o2.coupl)+1:length(o2.coupl)+length(o2.sigma)), tr(1:length(o2.coupl)), tr(length(o2.coupl)+length(o2.sigma)+1:end));
        end

        function r = mldivide(o1, o2)
            if ~isa(o1, 'DLSS_Hess')
                error('Type error!');
            end
            if ~isa(o2, 'DLSS_est')
                error('Type error!');
            end
            Hess = [o1.coupl o1.sc' o1.cz; o1.sc o1.sigma o1.sz; o1.cz' o1.sz' o1.zeta];
            grad = [o2.coupl; o2.sigma; o2.zeta];
            tr = Hess\grad;
            r = DLSS_est(tr(length(o2.coupl)+1:length(o2.coupl)+length(o2.sigma)), tr(1:length(o2.coupl)), tr(length(o2.coupl)+length(o2.sigma)+1:end));
        end

    end

end
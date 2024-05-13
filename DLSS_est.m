classdef DLSS_est% < handle


    properties
        sigma
        coupl
        zeta
    end


    methods

        function obj = DLSS_est(sigma, coupl, zeta)
            arguments
                sigma = []
                coupl = []
                zeta = []
            end
            obj.sigma = sigma;
            obj.coupl = coupl;
            obj.zeta = zeta;
        end

        function r = plus(o1, o2)
            if isa(o1, 'DLSS_est') && isa(o2, 'DLSS_est')
                s = o1.sigma + o2.sigma;
                c = o1.coupl + o2.coupl;
                z = o1.zeta + o2.zeta;
                r = DLSS_est(s, c, z);
            elseif isa(o1, 'DLSS_est')
                s = o1.sigma + o2;
                c = o1.coupl + o2;
                z = o1.zeta + o1;
                r = DLSS_est(s, c, z);
            else
                s = o2.sigma + o1;
                c = o2.coupl + o1;
                z = o2.zeta + o1;
                r = DLSS_est(s, c, z);
            end
        end

        function r = minus(o1, o2)
            if isa(o1, 'DLSS_est') && isa(o2, 'DLSS_est')
                s = o1.sigma - o2.sigma;
                c = o1.coupl - o2.coupl;
                z = o1.zeta - o2.zeta;
                r = DLSS_est(s, c, z);
            elseif isa(o1, 'DLSS_est')
                s = o1.sigma - o2;
                c = o1.coupl - o2;
                z = o1.zeta - o2;
                r = DLSS_est(s, c, z);
            else
                s = o1 - o2.sigma;
                c = o1 - o2.coupl;
                z = o1 - o2.zeta;
                r = DLSS_est(s, c, z);
            end
        end

        function r = uminus(o1)
            s = -o1.sigma;
            c = -o1.coupl;
            z = -o1.zeta;
            r = DLSS_est(s, c, z);
        end

        function r = mtimes(o1, o2)
            if isa(o1, 'DLSS_est')
                s = o1.sigma*o2;
                c = o1.coupl*o2;
                z = o1.zeta*o2;
                r = DLSS_est(s, c, z);
            elseif isa(o2, 'DLSS_est')
                s = o1*o2.sigma;
                c = o1*o2.coupl;
                z = o1*o2.zeta;
                r = DLSS_est(s, c, z);
            end
        end

        function r = norm(o1)
            r = sqrt(sum(o1.sigma.^2) + sum(o1.coupl.^2) + sum(o1.zeta.^2));
        end

    end

end
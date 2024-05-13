function p = inters(g1, g2, h1, h2)

    v1 = g2-g1;
    v2 = h2-h1;
    
    detm = (v1(2)*v2(1)-v1(1)*v2(2));
    if abs(detm) < 1e-10%the lines are parallel
        t = (h2-g1)./v1;
        if abs(t(1)-t(2)) < 1e-6%h2 is on the same line
            if t(1) >= 0 && t(1) <= 1%is h2 on the line segment?
                p = h2;
            else
                t = (h1-g1)./v1;
                t2= (g1-h1)./v2;
                if t(1) >= 0 && t(1) <= 1%is h1 on the line segment?
                    p = h1;
                elseif t2(1) >= 0 && t2(1) <= 1
                    p = g1;
                else
                    p = nan;
                end
            end
        else
            p = nan;
        end
        return;
    end
    t = (v2(2)*(g1(1)-h1(1)) + v2(1)*(h1(2)-g1(2)))/detm;

    if t < 0 || t > 1
        p = nan;
    else
        p = g1 + t*v1;
    end

end
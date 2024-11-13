function [Ln, dLn, ddLn] = LagrangePolynomial2(v, x, y)
    M = length(x);
    kmin = -2; kmax = 2;
    Ln(length(y)) = 0; 
    dLn(length(y)) = 0;
    ddLn(length(y)) = 0;
    
    for k = 1:length(y)
        q = -1;  
        min_dist = 2*(x(M) - x(1));

        for i=1:M
            if abs(y(k) -  x(i)) < min_dist
               min_dist = abs(y(k) -  x(i));
                q = i;
            end
        end
        idx = (q + kmin:q + kmax);
        if idx(1) < 1
            delta1 =  1 - idx(1);
            idx = idx + delta1;
        end
        if idx(end) > M
            delta2 = M - idx(end);
            idx = idx +delta2;
        end

        x1 = x(idx);
        for i=1:length(x1)
            p = 1;
            for j=1:length(x1)
                if j~=i
                    p = p * (y(k) - x1(j))/(x1(i) - x1(j));
                end       
            end
            Ln(k)= Ln(k) + p*v(idx(i));
        end

        for i = 1:length(x1)
            p = 1;
            for j = 1:length(x1)
                 if j~=i
                    p = p * (x1(i) - x1(j));
                end   
            end
            p = v(idx(i)) / p;
            sum = 0;
            for j = 1:length(x1)
                if j~=i
                    p2 = 1;
                    for l = 1:length(x1)
                        if l~=i && l~=j
                            p2 = p2 * (y(k) - x1(l)); 
                        end
                    end
                    sum = sum + p2;
                end
            end
            dLn(k) = dLn(k) + p * sum;
        end

        for i = 1:length(x1)
            p = 1;
            for j = 1:length(x1)
                 if j~=i
                    p = p * (x1(i) - x1(j));
                end   
            end
            p = v(idx(i)) / p;
            sum = 0;
            for j = 1:length(x1)
                if j~=i
                    sum2 = 0;
                    for l = 1:length(x1)
                        if l~=i && l~=j
                            p2 = 1;
                            for m = 1:length(x1)
                                if m~=i && m~=j && m~=l
                                    p2 = p2 * (y(k) - x1(m));
                                end
                            end
                            sum2 = sum2 + p2;
                        end
                    end
                    sum = sum + sum2;
                end
            end
            ddLn(k) = ddLn(k) + p * sum;
        end 
    end 
end 
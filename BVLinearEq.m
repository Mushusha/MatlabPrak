function v = BVLinearEq(y_n, ddy_n, h, M, cL, cR, cL_n, dcL_n, dcR_n)    
    for m = M:-1:2
        a(m) = 1/h^2;
        b(m) = -2/h^2 - 12 * y_n(m) + 4;
        c(m) = 1/h^2;
        d(m) = -ddy_n(m) - 4*y_n(m) + 6*y_n(m)*y_n(m);

    end
    a(1) = 0; b(1) = 1/2 + cL/h; c(1) = 1/2 - cL/h; d(1) = cL * dcL_n - cL_n;
    a(M+1) = -1/h; b(M+1) = 1/h; c(M+1) = 0; d(M+1) = cR - dcR_n;
    v = Progonka(a, b, c, d);
end
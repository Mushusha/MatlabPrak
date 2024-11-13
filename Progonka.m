function xn = Progonka(a, b, c, d)
    yn = 0 * a;
    yn(1) = b(1);
    A = 0*a;
    B = 0*a;
    A(1) = -c(1)/yn(1); B(1) = d(1)/yn(1);
    [l,n] = size(A);
    for m = 2:n-1
        yn(m) = b(m) + a(m)*A(m-1);
        A(m) = -c(m)/yn(m); B(m) = (d(m) - a(m)*B(m-1)) / yn(m);
    end
    yn(n) = b(n) + a(n)*A(n-1);
    B(n) = (d(n) - a(n)*B(n-1)) / yn(n);

    xn(n) = B(n);
    for m = (n-1):-1:1  
        xn(m) = A(m)*xn(m+1) + B(m);
    end
end
clear
format long
syms z

xL = 0.1; 
xR = 3 * pi / 4; 
eps = 1e-4; 
alpha = 0;
cL = 0.5 * tan(xL); 
cR = 2 * cot (xR) * func(xR);
iter_array = ['iter01','iter02','iter03','iter04','iter05','iter06','iter07','iter08','iter09','iter10','iter11','iter12','iter13','iter14','iter15'];
mesh_array = ['mesh01','mesh02','mesh03','mesh04','mesh05','mesh06','mesh07','mesh08','mesh09','mesh10','mesh11','mesh12','mesh13','mesh14','mesh15','mesh16','mesh17','mesh18'];

v2h = 1;
x2h = 1;

M = 2^10;
h = (xR - xL) / M;
x_n = xL - h / 2 + (0:M + 1) * h;
y_n = (1 - alpha) * func(x_n) + alpha * (cR * x_n - cR * xL + cR * cL);
figure(1);
hold off
plot(x_n, y_n);
hold on
legend('y_{ini}');
plot(x_n, 1 ./ func(x_n),'DisplayName','y_a');
figure(2);
hold off
title(['Convergence of newtotian iterations']);
hold on
figure(3);
hold off
legend();
title(['NewtonianItertion=',num2str(1)]);
hold on

for steps = 1:15
    fprintf('Step %d==============================\n', steps)
    delta2h = 2;
    delta1h = 2;
    for mesh = 1:16
        M = 2^(mesh+4);
        h = (xR - xL) / M;
        x = xL - h / 2 + (0:M + 1) * h;
        [Ly_n, dLy_n, ddLy_n] = LagrangePolynomial2(y_n, x_n, x);
        [cL_n,dcL_n,tmp2] = LagrangePolynomial2(y_n, x_n, xL);
        [tmp3,dcR_n,tmp4] = LagrangePolynomial2(y_n, x_n, xR);

        v1h = BVLinearEq(Ly_n, ddLy_n, h, M+1, cL, cR, cL_n, dcL_n, dcR_n);
        if mesh>1
            L_v2 = LagrangePolynomial2(v2h, x2h, x);
            delta1h = max(abs(v1h - L_v2));
        end
        ratio = delta2h/delta1h;
        v2h = v1h;
        x2h = x;
        delta2h = delta1h;
        fprintf('M = %d, ratio = %f, delta = %f\n', M, ratio, delta1h);
        if delta1h < eps/2
            fprintf('Mesh = %d\n', mesh);
            break
        end
        figure(3);
        legend();
        if(mesh~=1)
            hold on
        end
        title(['NewtonianIteration= ',num2str(steps),' mesh= ',num2str(mesh),' delta= ',num2str(delta1h),' ratio= ',num2str(ratio)]);
        plot(x, v1h,'DisplayName',mesh_array(1 + (mesh-1)*6:6 + (mesh-1)*6));
        pause(0.5);
    end

    if delta1h > eps/2
        break
    end
    y_n = y_n + LagrangePolynomial2(v2h, x2h, x_n);
    figure(1);
    title(['NewtonianIter = ',num2str(steps),' ||v|| = ',num2str(max(abs(v1h)))]);
    plot(x_n, y_n,'.','DisplayName',iter_array(1 + (steps-1)*6:6 + (steps-1)*6));
    figure(2);
    xlim([0, 10]);
    plot(steps,log10(max(abs(v1h))),'o');
    pause(0.2);
    fprintf('||v|| = %f\n', max(abs(v1h)));
    if(max(abs(v1h)) < eps/2)
        fprintf('Steps = %d\n', steps);
        break
    end
    figure(3)
    hold off

end
fprintf('||y_n - y_a|| = %f\n', max(abs(y_n - func(x_n))));
if(max(abs(y_n - func(x_n)) < eps))
    fprintf('Converged for alpha = %f\n', alpha);
else
    fprintf('Didnt converge for alpha = %f\n', alpha);
end
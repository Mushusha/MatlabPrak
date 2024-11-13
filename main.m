clear
format long

xL = 0.1; 
xR = 3 * pi / 4; 
eps = 1e-4; 
alpha = 0.49; %alpha [-0.17; 0.49]
fprintf('alpha = %e\n', alpha);
cL = -0.5 * tan(xL); 
cR = -2 * cot (xR) * func(xR);
maxsteps = 15;
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
plot(x_n, y_n,'g');
ymin=min(y_n);
ymax=max(y_n);

hold on
legend('y_{ini}');
yanalit=func(x_n);
plot(x_n, yanalit,'r','DisplayName','y_a');
ymin=min(ymin,min(yanalit))-1;
ymax=max(ymax,max(yanalit))+1;
axis([xL xR ymin, ymax]);
figure(2);
hold off
plot([xL xR],[0 0],'w')
title(['Correction v: NewtonianItertion=',num2str(1)]);
figure(3);
hold off
plot([0 maxsteps],[log10(eps) log10(eps)],'k--')
hold on
plot([1 1],[2 log10(eps)-1],'w')
title(['Convergence of newtotian iterations']);
legend();

for steps = 1:maxsteps
    fprintf('Step %d\n', steps)
    delta2h = 2;
    delta1h = 2;
    figure(2)
    hold off
    for mesh = 1:11
        tstart=cputime;
        M = 2^(mesh+4);
        h = (xR - xL) / M;
        x = xL - h / 2 + (0:M + 1) * h;
        [Ly_n, dLy_n, ddLy_n] = LagrangePolynomial2(y_n, x_n, x);
        [cL_n,dcL_n,tmp2] = LagrangePolynomial2(y_n, x_n, xL);
        [tmp1,dcR_n,tmp4] = LagrangePolynomial2(y_n, x_n, xR);

        v1h = BVLinearEq(Ly_n, ddLy_n, h, M+1, cL, cR, cL_n, dcL_n, dcR_n);
        if mesh>1
            L_v2 = LagrangePolynomial2(v2h, x2h, x);
            delta1h = max(abs(v1h - L_v2));
        end
        ratio = delta2h/delta1h;
        v2h = v1h;
        x2h = x;
        delta2h = delta1h;
        fprintf('mesh =%2d M = %6d, ratio = %f, delta = %e t = %g sec\n',...
            mesh, M, ratio, delta1h,cputime-tstart);
        if delta1h < eps/2
            fprintf('Mesh = %d\n', mesh);
            break
        end
        figure(3);
        legend();
        if(mesh~=1)
            hold on
        end
        figure(2)
        plot(x, v1h,'DisplayName',mesh_array(1 + (mesh-1)*6:6 + (mesh-1)*6));
        xlim([xL, xR]);
        title(['Correction v: NewtonianIteration = ',num2str(steps),' mesh = ',num2str(mesh),' delta = ',num2str(delta1h),' ratio = ',num2str(ratio)]);
        hold on
        pause(0.5);
    end 

    if delta1h > eps/2
        fprintf('Didn''t converge for steps = %d: delta1h = %e > eps/2 = %e\n',...
            steps, delta1h, eps/2);
    end
    y_n = y_n + LagrangePolynomial2(v2h, x2h, x_n);
    figure(1);
    if fix(steps/2)*2 == steps
        color='c';
    else
        color='m';
    end
    plot(x_n, y_n,color);
    axis([xL xR ymin, ymax]);
    title(['NewtonianIter = ',num2str(steps),' ||v|| = ',num2str(max(abs(v1h)))]);
    figure(3)
    plot(steps,log10(max(abs(v1h))),'.','MarkerSize',20);
    title(['Convergence of newtotian iterations with \alpha =',num2str(alpha)]);
    pause(0.2);
    fprintf('alpha = %e, ||v|| = %e\n', alpha, max(abs(v1h)));
    if(max(abs(v1h)) < eps/2)
        fprintf('Steps = %d\n', steps);
        break
    end
    pause(1)

end
fprintf('||y_n - y_a|| = %f\n', max(abs(y_n - func(x_n))));
if(max(abs(y_n - func(x_n))) < eps)
    fprintf('Converged for alpha = %f\n', alpha);
else
    fprintf('Didn''t converge for alpha = %f\n', alpha);
end
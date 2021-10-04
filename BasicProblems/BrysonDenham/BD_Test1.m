toms t tf
p = tomPhase('p', t, 0, tf, 30); setPhase(p);
tomStates x1 x2
x1max = 1/9; x0 = {tf == 0.5};

constr = {0.001 <= tf <= 50
    collocate({0 <= x1 <= x1max; -10 <= x2 <= 10})
    initial({x1 == 0; x2 == 1}); final({x1 == 0; x2 == -1})
    collocate(dot(x1) == x2)};

options = struct;
options.name = 'Bryson Denham Short';
solution = ezsolve(integrate(0.5*dot(x2).^2), constr, x0, options);
t  = subs(collocate(t),solution);
x1 = subs(collocate(x1),solution); x2 = subs(collocate(x2),solution);
figure(1)
plot(t,x1,'*-',t,x2,'*-');
legend('x1','x2');
title('Bryson Denham Short state variables');
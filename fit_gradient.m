function [t, fval] = fit_gradient(particle, template, lmobj)
ds = lmobj.S(2)-lmobj.S(1);
sigmaedges = [lmobj.S - ds/2, inf];
options = optimoptions(@fminunc,'Display','none');
particle(:,3) = discretize(particle(:,3),sigmaedges);
func = @(p) costfunc(p, particle, template, lmobj.X, lmobj.Y, 1:size(template,3));
[p, fval] = fminunc(func, [0,0,0], options);
t = affine2d([cos(p(3)) sin(p(3)) 0; -sin(p(3)) cos(p(3)) 0; p(1) p(2) 1]);

function ll = costfunc(p, data, V, X, Y, Z)
dt = p(3);
t = affine2d([cos(dt) sin(dt) 0; -sin(dt) cos(dt) 0; p(1) p(2) 1]);
d2 = t.transformPointsForward(data(:,1:2));
ll = - sum(interp3(X, Y, Z, V, d2(:,1), d2(:,2), data(:,3), 'linear', log(eps(0))));

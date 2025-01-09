function [XX,CC] = trefoilpatch(r1,r2,r3,d1,d2);

%     r1 = 1.0
%     r2 = 3.0
%
%     Bx [u_, v_] := (r2 + r1 * Cos[u/2.0]) * Cos[u/3.0]
%     By [u_, v_] := (r2 + r1 * Cos[u/2.0]) * Sin[u/3.0]
%     Bz [u_, v_] := r1 * Sin[u/2.0]
%
%     x [u_, v_] := N[Bx [u, v]] + r1 * Cos[u/3.0] * Cos[v - Pi]
%     y [u_, v_] := N[By [u, v]] + r1 * Sin[u/3.0] * Cos[v - Pi]
%     z [u_, v_] := N[Bz [u, v]] + r1 * Sin[v - Pi]
%
%     ParametricPlot3D[N[{x[u, v], y[u, v], z[u, v]}],
%  	 {u, 0, 12 Pi}, {v, 0, 2 Pi}, PlotPoints -> {144, 12}]

%     r1 = 1.0
%     r2 = 3.0

u = 0:d1:12*pi; 
v = 0:d2:2*pi;
[u v] = meshgrid(u,v);
u = u(:)';
v = v(:)';

Bx = (r3 + r2 * cos(u/2.0)) .* cos(u/3.0);
By = (r3 + r2 * cos(u/2.0)) .* sin(u/3.0);
Bz = r2 * sin(u/2.0);


x = Bx + r1 * cos(u/3.0) .* cos(v - pi);
y = By + r1 * sin(u/3.0) .* cos(v - pi);
z = Bz + r1 * sin(v - pi);

XX = [x;y;2*z];
CC = [u;v];
%     ParametricPlot3D[N[{x[u, v], y[u, v], z[u, v]}],
%  	 {u, 0, 12 Pi}, {v, 0, 2 Pi}, PlotPoints -> {144, 12}]

% PINEQUALITYSURF  help with proving p-monotone inequality in Appendix by plotting

f = @(t,s,p) (1 - (t.^(p-1) + t) .* s + t.^p) ./ (1 - 2*s.*t + t.^2).^(p/2);

tt = 0:.01:1;
ss = -1:.01:.99;
[ttt,sss]=meshgrid(tt,ss);

surf(sss,ttt,f(ttt,sss,2.5)), colorbar, shading('interp')


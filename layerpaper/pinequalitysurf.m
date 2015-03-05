function pinequalitysurf(p)
% PINEQUALITYSURF  help with proving p-monotone inequality in Appendix by plotting

f = @(t,s,p) (1 - (t.^(p-1) + t) .* s + t.^p) ./ (1 - 2*s.*t + t.^2).^(p/2);

tt = 0:.01:1;
ss = -1:.01:.99;
[ttt,sss]=meshgrid(tt,ss);

figure(1)
surf(sss,ttt,f(ttt,sss,p)), colorbar, shading('interp')
title('f(t,s)')

figure(2)
g = @(t,s,p) s .* (2-p) .* t .* (t.^(p-2) + 1) + ((p-1) * (t.^p+1) - (t.^(p-2) + t.^2));
surf(sss,ttt,g(ttt,sss,p)), colorbar, shading('interp')
title('g(t,s) = df/ds / (positive factor)')

figure(3)
plot(tt,(p-1)*(tt.^p+1),'r',tt,(p-2)*(tt.^(p-1)+tt)+tt.^(p-2)+tt.^2,'b')
legend('should be above','should be below','location','southeast')


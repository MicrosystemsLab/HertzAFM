function sse = funHertz(x)
global R xdata ydata F
F = (4/3)*x(1)*xdata.^(3/2)*sqrt(R);
error = ydata-F;
sse = sum(error.^2); 
% from Radmacher 2007 and Carl 2008
%x is Ereduced

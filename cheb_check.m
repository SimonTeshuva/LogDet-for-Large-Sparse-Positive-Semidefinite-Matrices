close all;
clear all;
clc;

delta = 1e-3;
n = 15;
x = chebfun('x');
p = log((1-((1-2*delta).*x+1)/2)); % using this d we get a complex chebyshev fn
c = chebpoly(p, n); % this is not the funcion you are looking for
c = fliplr(c); 

c = fliplr(c); 
x2 = -1:0.01:1;
y = log((1-((1-2*delta).*x2+1)/2));
clf, plot(c)
hold on
plot(x2,y)
hold off
legend('approx', 'exact')
figure (2)
plotcoeffs(p)

cp = zeros(size(x2));
for counter = 1:(length(c)-1)
    cp = cp + c(counter)*x2.^(counter-1)
end
cp = cp+c(end)
figure (3)
plot(cp)
plot(x2, cp)
plot(x2, cp)
hold on 
plot(x2,y)
hold off
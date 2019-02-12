function myfunc(x)
    x1, x2 = x[1], x[2]
    return x1^4 - 2*x2*x1^2 + x2^2 + x1^2 - 2*x1 + 5
end

lb = [-3.0, -1.0]
ub = [2.0, 6.0];

particles = Particles(100, lb, ub)
pso(particles, myfunc, minstep=1e-3)
# ou posso fazer a função criar um e depor retornar
particles2 = pso(myfunc, lb, ub, minstep=1e-3);

mínimo em [1,1] com valor 4
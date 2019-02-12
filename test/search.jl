# melhorar esse teste aqui

particles = Particles(5, zeros(3), ones(3).*5)
particles.x = [fill(a,3) for a in 1:5]
particles.fx = [ a*0.1 for a in 5:-1:1]
particles.p = copy(particles.x)
particles.fp = copy(particles.fx) .* 0.1
particles.g = copy(particles.x)
particles.fg = copy(particles.fx)

findBestLocal(particles, 1)

[0.01, 0.03, 0.02, 0.01, 0.01]

[555
333
444
555
555]

findBestGlobal(particles)
particles.fg =  0.010000000000000002
0.010000000000000002
0.010000000000000002
0.010000000000000002
0.010000000000000002

particles.g = [5.0, 5.0, 5.0]
[5.0, 5.0, 5.0]
[5.0, 5.0, 5.0]
[5.0, 5.0, 5.0]
[5.0, 5.0, 5.0]
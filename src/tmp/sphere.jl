    # known locations
    front = (0.0, 1.0, 0.0)
    right = (1.0, 0.0, 0.0)
    back = (0.0, -1.0, 0.0)
    left = (-1.0, 0.0, 0.0)
    top = (0.0, 0.0, 1.0)

# sphere
n = 20;
if nargs > 0, n = args{1} ;end

% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;

x = cosphi*cos(theta);
y = cosphi*sintheta;
z = sin(phi)*ones(1,n+1);


    # draw spherical head
    resolution = 100
    u, v = np.mgrid[0 : 2 * np.pi : resolution, 0 : np.pi : resolution]

grid = Iterators.product(
                         u = collect(range(0, 2 * pi, length=resolution))
                         v = collect(range(0, pi, length=resolution))
                         );
collect.(grid)

n = 20

theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos.(phi); cosphi[1] = 0; cosphi[n+1] = 0;
sintheta = sin.(theta); sintheta[1] = 0; sintheta[n+1] = 0;

x = @. cosphi*cos(theta);
y = @. cosphi*sintheta;
z = @. sin(phi) * ones(1, n+1);

using GLMakie

x, y = collect(-8:0.5:8), collect(-8:0.5:8)
z = [sinc(√(X^2 + Y^2) / π) for X ∈ x, Y ∈ y]

wireframe(x, y, z, axis=(type=Axis3,), color=:black)

N = 32
u = range(0, stop=2π, length=N)
v = range(0, stop=π, length=N)
x = cos.(u) .* sin.(v)'
y = sin.(u) .* sin.(v)'
z = repeat(cos.(v)',outer=[N, 1])

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

wireframe  (xs, ys, zs, axis=(type=Axis3,))

using GLMakie

fig = Figure()
ax = Axis3(fig[1, 1]; aspect=(1, 1, 1), perspectiveness=0.5)
surface(ax, xs, ys, zs)

for idx in 1:resolution
    GLMakie.scatter!(ax, x[idx], y[idx], z[idx])
end
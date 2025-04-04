export plot_connectivity_circle

m = rand(10, 10)
m = rand(10, 10)
labels = string.(collect('A':'Z')[1:size(m, 1)])
t = linspace(-pi, pi, size(m, 1))
pos_x = zeros(size(m, 1))
pos_y = zeros(size(m, 1))
pos_col = zeros(size(m, 1))
for idx in axes(m, 1)
  pos_x[idx] = cos(t[idx])
  pos_y[idx] = sin(t[idx])
  pos_col[idx] = idx
end
using Plots

p = Plots.plot(aspect_ratio = :equal,
               framestyle=:none,
               legend=false)
for idx in axes(m, 1)
  p = Plots.scatter!((pos_x[idx], pos_y[idx]))
end

P = (pos_x[1], pos_y[1])
Q = (pos_x[4], pos_y[4])

m_PQ = ((P[1] + Q[1]) / 2, (P[2] + Q[2]) / 2)
slope_PQ = (Q[2] - P[2]) / (Q[1] - P[1])
slope_bisector = - 1/slope_PQ

alpha = 1 / (P[1]^2 + P[2]^2)
beta = 1 / (Q[1]^2 + Q[2]^2)
P2 = (alpha * P[1], alpha * P[2])
Q2 = (beta * Q[1], beta * Q[2])

# Draw connections on the Poincare hyperbolic disk.
#
# Equation of the circles on the disk:
# x^2 + y^2 
# + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x 
# - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
# where u and v are points on the boundary.
#
# Standard form of equation of a circle
# (x - x0)^2 + (y - y0)^2 = r^2
#
# Therefore we can identify
# x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
# y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
# r^2 = x0^2 + y0^2 - 1

u = zeros(size(m, 1))
v = zeros(size(m, 1))
for idx in axes(m, 1)
  u  = [cos()tpos_x = zeros(size(m, 1))
pos_y = zeros(size(m, 1)), sin(t(row(i)))
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0 
              % ensure the arc is within the unit disk
              theta = [linspace(max(thetaLim),pi,50),...
                       linspace(-pi,min(thetaLim),50)].';
            else
              theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            
            this.Node(row(i)).Connection(end+1) = line(...
              r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', lineWidth(i),...
              'Color', this.ColorMap(row(i),:),...
              'PickableParts','none');


arc = Plots.partialcircle(1, 1//2 * pi)
for idx in 1:length(arc)
  arc[idx] = (arc[idx][1] + 1, arc[idx][2])
end
Plots.plot!(arc)
Plots.plot(p)
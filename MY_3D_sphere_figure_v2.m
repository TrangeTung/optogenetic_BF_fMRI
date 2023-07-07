function F = MY_3D_sphere_figure_v2(V,ColorVec)


F = figure('Position', [27 72 1071 906]);
t = linspace(0,pi,50);
p= linspace(0,2*pi,50);
[theta,phi] = meshgrid(t,p);
x = sin(theta).*sin(phi);
y = sin(theta).*cos(phi);
z = cos(theta);
h=surf(x,y,z,'facecolor','k','linestyle','none');
hold on
set(h,'FaceColor','interp');
set(h,'EdgeColor','none');
axis equal;
zlabel('NMF C3');
xlabel('NMF C1');
ylabel('NMF C2');
MyMap = gray(125);
% MyMap = [ones(125,1),(1:125)'/125,(1:125)'/125];
colormap(MyMap);
light;
set(h,'FaceLighting','gouraud');
set(h,'DiffuseStrength',.3);
set(h,'AmbientStrength',.5);
set(h,'SpecularStrength',.4);
set(h,'SpecularExponent',2);
plot3([0, 1.5] ,[0 0], [0 0], 'k', 'linewidth',1.5)
plot3([0 0],[0, 1.5] , [0 0], 'k', 'linewidth',1.5)
plot3([0 0], [0 0],[0, 1.5] , 'k', 'linewidth',1.5)
axis equal; grid off
alpha(0.8);
light('position',[10 30 30],'color','r')
view([135 45]);

theta0 = 0;
theta1 = theta0 + 2*pi;
granularity = 500;
thetas = linspace(theta0,theta1,ceil((theta1-theta0)/(2*pi) * (granularity+1)));
[X,Y] = pol2cart(thetas,1);
coord = [X; Y];
sd1 = 1;
sd2 = 1;
coord = diag([sd1 sd2]) * coord;
ang = 0;
coord = [cos(-ang) sin(-ang); -sin(-ang) cos(-ang)]*coord;
plot3(coord(1,:),coord(2,:),zeros(size(coord(1,:))),'w-','linewidth',1.5);
plot3(coord(1,:),zeros(size(coord(1,:))),coord(2,:),'w--','linewidth',1.5);
plot3(zeros(size(coord(1,:))),coord(1,:),coord(2,:),'w--','linewidth',1.5);

for vl=1:numel(V)
   Vec = V{vl};
   V0 = Vec([1,3,4])./norm(Vec([[1,3,4]]),2)*1.2;
   plot3([0 V0(1)],[0 V0(2)],[0 V0(3)],'-o','linewidth',2,'color',ColorVec{vl}/255);
   hold on
end

Vm = mean(cell2mat(V),2)*1.2;
V0 = Vm([1,3,4])./norm(Vm([[1,3,4]]),2)*2.5;
plot3([0 Vm(1)],[0 V0(2)],[0 V0(3)],'linewidth',4,'color','k');
end
clear all; close all; clc;

rho_t = 0.5;

Xe= [25.3409 22.05 24.375 28.53158 25.3409 22.05 24.375 28.53158; -36.29172 -38.19172 -42.21874 -40.80777 -36.29172 -38.19172 -42.21874 -40.80777; 9.166667 9.166667 9.166667 9.166667 13.33333 13.33333 13.33333 13.33333]
rho_e = [0.375, 0.5, 1.0, 0.8333333333333334, 0.375, 0.5, 1.0, 0.8333333333333334]
x = [15.0, -60.0, -19.5]

plot3(Xe(1,:), Xe(2,:), Xe(3,:),'ro'); grid on; hold on;
plot3(x(1), x(2), x(3),'gx'); grid on;

xi = -1.0:0.1:1.0;

[Xi1, Xi2, Xi3] = meshgrid(xi,xi,xi);

N = sfce(Xi1, Xi2, Xi3);

X = zeros(size(Xi1));
Y = zeros(size(Xi1));
Z = zeros(size(Xi1));
rho = zeros(size(Xi1));

for i = 1:8
    X = X + N(:,:,:,i) * Xe(1,i);
    Y = Y + N(:,:,:,i) * Xe(2,i);
    Z = Z + N(:,:,:,i) * Xe(3,i);
    rho = rho + N(:,:,:,i) * rho_e(i);
end

figure
isosurface(Xi1, Xi2, Xi3,rho, 0.5)
figure
isosurface(X, Y, Z,rho, 0.5)

function N = sfce(xi1, xi2, xi3)    
    N1 = -1/8*(xi1-1).*(xi2-1).*(xi3-1);
    N2 =  1/8*(xi1+1).*(xi2-1).*(xi3-1);
    N3 = -1/8*(xi1+1).*(xi2+1).*(xi3-1);
    N4 =  1/8*(xi1-1).*(xi2+1).*(xi3-1);

    N5 =  1/8*(xi1-1).*(xi2-1).*(xi3+1); %!!!!
    N6 = -1/8*(xi1+1).*(xi2-1).*(xi3+1); %!!!!
    N7 =  1/8*(xi1+1).*(xi2+1).*(xi3+1);
    N8 = -1/8*(xi1-1).*(xi2+1).*(xi3+1);
    N = cat(4, N1, N2, N3, N4, N5, N6, N7, N8);
end
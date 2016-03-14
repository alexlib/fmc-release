function VORTICITY_GRAD = calculateVorticity_socdiff(X, Y, U, V)

% Spacing between grid points
dx = X(1,2) - X(1,1);
dy = Y(2, 1) - Y(1, 1);

% Partial derivatives of velocity
dvdx = socdiff(V, dx, 2);
dudy = socdiff(U, dy, 1);

VORTICITY_GRAD = -1 * (dvdx - dudy);


end


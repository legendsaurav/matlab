function cmap = plasmaLike(n)
%PLASMALIKE  Small dependency-free plasma-style colormap sampler.
%   cmap = PLASMALIKE(n) returns an n x 3 RGB matrix.
    base = [13 8 135; 126 3 168; 204 71 120; 248 149 64; 240 249 33] / 255;
    x  = linspace(0, 1, size(base,1));
    xi = linspace(0, 1, n);
    cmap = [interp1(x, base(:,1), xi, 'linear').', ...
            interp1(x, base(:,2), xi, 'linear').', ...
            interp1(x, base(:,3), xi, 'linear').'];
end

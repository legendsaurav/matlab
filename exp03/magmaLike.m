function cmap = magmaLike(n)
%MAGMALIKE  Small dependency-free magma-style colormap sampler.
%   cmap = MAGMALIKE(n) returns an n x 3 RGB matrix.
    base = [0 0 4; 81 18 124; 183 55 121; 252 137 97; 252 253 191] / 255;
    x  = linspace(0, 1, size(base,1));
    xi = linspace(0, 1, n);
    cmap = [interp1(x, base(:,1), xi, 'linear').', ...
            interp1(x, base(:,2), xi, 'linear').', ...
            interp1(x, base(:,3), xi, 'linear').'];
end

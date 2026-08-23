function cmap = viridisLike(n)
%VIRIDISLIKE  Small dependency-free viridis-style colormap sampler.
%   cmap = VIRIDISLIKE(n) returns an n x 3 RGB matrix. Used instead of
%   MATLAB's built-in viridis() purely so this code also runs unmodified
%   on older MATLAB releases / Octave installs that lack it.
    base = [68 1 84; 59 82 139; 33 145 140; 94 201 98; 253 231 37] / 255;
    x  = linspace(0, 1, size(base,1));
    xi = linspace(0, 1, n);
    cmap = [interp1(x, base(:,1), xi, 'linear').', ...
            interp1(x, base(:,2), xi, 'linear').', ...
            interp1(x, base(:,3), xi, 'linear').'];
end

function Bout = mlpcr2_bootfun(X, Y)

% takes subjXblockXvoxel data as X
% takes xubjXblock data as Y
nbetween = size(X, 1);
nwithin = size(X, 2);
nvox = size(X, 3);

Ynew = reshape(Y, [nbetween*nwithin, 1]);
Xnew = reshape(X, [nbetween*nwithin, nvox]);
subjIDs = repelem(1:nbetween, nwithin)';


[B, Bb, Bw, pc_b, sc_b, pc_w, sc_w] = mlpcr2(Xnew, Ynew, 'subjIDs', subjIDs);
Bout = [Bb' Bw'];

end
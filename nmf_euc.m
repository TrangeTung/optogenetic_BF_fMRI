function [W,H,errs,loss] = nmf_euc(V, r)


if min(min(V)) < 0
  error('Matrix entries can not be negative');
end
if min(sum(V,2)) == 0
  error('Not all entries in a row can be zero');
end

[m,n] = size(V);
W = rand(m,r);
H = rand(r,n);

niter = 1000;

myeps = 1e-10;

errs = zeros(niter,1);

for t = 1:niter

   W = W .* ( (V*H') ./ max(W*(H*H'), myeps) ); 
   %W = normalize_W(W,1);
   H = H .* ( (W'*V) ./ max((W'*W)*H, myeps) );

   loss = sum((V-W*H).^2);
   errs(t) = sum(sum(loss));
end
1;
end
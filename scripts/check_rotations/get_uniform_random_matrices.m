function M = get_uniform_random_matrices( N )
% M = get_uniform_random_matrices( N )

for i = 1:N
 [Q, R] = qr(randn(3));
 Q = Q*diag( sign(diag(R)));
 if det(Q) < 0; Q(:,1) = -Q(:,1); end;
 M(:,:,i) = Q;
end

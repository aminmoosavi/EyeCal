function E = rieszEnvelope(X)
% RIESZENVELOPE  Local envelope via the monogenic signal
%   E = rieszEnvelope(X) returns sqrt( X.^2 + R1.^2 + R2.^2 ),
%   where [R1,R2] are the Riesz transform of X.
% compute Riesz components
[M,N] = size(X);
F = fft2(X);
ux = ifftshift((-floor(N/2):ceil(N/2)-1)/N);
uy = ifftshift((-floor(M/2):ceil(M/2)-1)/M).';
[U,V] = meshgrid(ux,uy);
K = sqrt(U.^2+V.^2)+eps;
H1 = -1i*(U./K);
H2 = -1i*(V./K);
R1 = ifft2(F.*H1);
R2 = ifft2(F.*H2);
% envelope = magnitude of [X, R1, R2]
E = sqrt( abs(X).^2 + abs(R1).^2 + abs(R2).^2 );
end
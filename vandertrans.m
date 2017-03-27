A=rand(4,4);
AA=A'*A
eig_AA=eig(AA)
vander_AA=vander(eig_AA)
T=flipud(vander_AA')
pinv(T)*AA*T
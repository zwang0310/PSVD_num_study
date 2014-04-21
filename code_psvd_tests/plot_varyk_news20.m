result_rsvdp1 = [400 5.211000e+01 5.206120e-01 1.292258e-02 
  800 1.606500e+02 4.485109e-01 8.878910e-03 
 1200 3.664400e+02 3.982745e-01 7.126325e-03 
 1600 6.169400e+02 3.583585e-01 6.024544e-03 
 2000 9.518200e+02 3.246936e-01 5.226722e-03];

result_rsvdp2 = [400 1.672300e+02 5.146084e-01 6.385328e-03 
  800 6.165000e+02 4.413467e-01 3.860387e-03 
 1200 1.431970e+03 3.906492e-01 2.638701e-03 
 1600 2.607420e+03 3.504370e-01 1.888187e-03 
 2000 4.249560e+03 3.168071e-01 1.349142e-03];

result_rsvdk = [400 3.406000e+02 5.155687e-01 1.064938e-02 
  800 1.360490e+03 4.422424e-01 6.835850e-03 
 1200 3.127800e+03 3.913763e-01 4.968223e-03 
 1600 6.188150e+03 3.509688e-01 3.697800e-03 
 2000 1.084534e+04 3.171628e-01 2.699751e-03];

result_propack = [400 2.956200e+02 5.140454e-01 5.237374e-12 
  800 1.538140e+03 4.408694e-01 5.675619e-12 
 1200 3.306630e+03 3.903150e-01 1.182718e-12 
 1600 3.787810e+03 3.502230e-01 1.008600e-12 
 2000 7.216120e+03 3.166776e-01 9.943501e-13];

result_svds = [400 8.571800e+02 5.140454e-01 1.736272e-14 
  800 4.584490e+03 4.408694e-01 1.617550e-13 
 1200 1.000578e+04 3.903150e-01 1.147148e-10 
 1600 1.919582e+04 3.502230e-01 1.569496e-14 
 2000 2.475836e+04 3.166776e-01 1.583645e-14];

result_lmsvd = [400 1.923110e+03 5.140454e-01 1.033993e-06 
  800 6.553330e+03 4.408694e-01 1.414608e-06 
 1200 1.612318e+04 3.903150e-01 1.791815e-06 
 1600 2.436323e+04 3.502230e-01 2.091676e-06 
 2000 3.460495e+04 3.166776e-01 2.253753e-06];

result_bchdav = [400 2.437900e+02 5.140454e-01 4.836097e-17
    800 6.649400e+02 4.408694e-01 4.232492e-17
    1200 1.305200e+03 3.903150e-01 3.848468e-17
    1600 2.701800e+03 3.502230e-01 4.273895e-17
    2000 4.258450e+03 3.166776e-01 4.664119e-17];

kvals = result_rsvdp1(:,1);

t_rsvdp1 = result_rsvdp1(:,2);
t_rsvdp2 = result_rsvdp2(:,2);
t_rsvdk  = result_rsvdk(:,2);
t_propack = result_propack(:,2);
t_svds = result_svds(:,2);
t_lmsvd = result_lmsvd(:,2);
t_bchdav = result_bchdav(:,2);

merr_rsvdp1 = result_rsvdp1(:,3);
merr_rsvdp2 = result_rsvdp2(:,3);
merr_rsvdk  = result_rsvdk(:,3);
merr_propack = result_propack(:,3);
merr_svds = result_svds(:,3);
merr_lmsvd = result_lmsvd(:,3);
merr_bchdav = result_bchdav(:,3);

verr_rsvdp1 = result_rsvdp1(:,4);
verr_rsvdp2 = result_rsvdp2(:,4);
verr_rsvdk  = result_rsvdk(:,4);
verr_propack = result_propack(:,4);
verr_svds = result_svds(:,4);
verr_lmsvd = result_lmsvd(:,4);
verr_bchdav = result_bchdav(:,4);

figure(1)
subplot(1,2,1)
plot(kvals,t_rsvdp2,'k-+',kvals,t_rsvdk,'b-+',kvals,t_propack,'m-o',...
    kvals,t_svds,'g-x',kvals,t_lmsvd,'c-*',kvals,t_bchdav,'r->')
legend('rsvd\_power (p=k q=2)','rsvd\_krylov (p=2 q=2)','PROPACK',...
    'svds (MATLAB)', 'LMSVD', 'bchdav')
xlabel('k: number of singular triplets')
ylabel('CPU time')
title('CPU time w.r.t. number of singular triplets')

subplot(2,2,2)
plot(kvals,merr_rsvdp2,'k-+',kvals,merr_rsvdk,'b-+',kvals,merr_propack,'m-o', ...
    kvals,merr_svds,'g-x',kvals,merr_lmsvd,'c-*',kvals,merr_bchdav,'r->')
xlabel('k: number of singular triplets')
ylabel('err\_mat')
title('err\_mat w.r.t. number of singular triplets')

subplot(2,2,4)
semilogy(kvals,verr_rsvdp2,'k-+',kvals,verr_rsvdk,'b-+',kvals,verr_propack,'m-o',...
    kvals,verr_svds,'g-x',kvals,verr_lmsvd,'c-*',kvals,verr_bchdav,'r->')
xlabel('k: number of singular triplets')
ylabel('max(err\_vec)')
title('max(err\_vec) w.r.t. number of singular triplets')
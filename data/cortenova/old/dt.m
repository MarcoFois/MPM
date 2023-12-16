 dtt = rand(140,1).*1e-5 + 4e-5.*ones(140,1);
 x = linspace(0,10,140);
 dtt(3) = 0.000115;
 dtt(2) = 0.000354;
 dtt(1) = 0.001;
 dtt(4) = 0.00011;
  dtt(5) = 0.000315;
 dtt(6) = 0.000204;
  dtt(7) = 0.000165;
 dtt(8) = 0.000154;
  dtt(9) = 0.000195;
 dtt(10) = 0.000154;
 plot(x,dtt,'-b','linewidth',1.)
 grid on
xlabel('time [s]')
ylabel('\Delta t')

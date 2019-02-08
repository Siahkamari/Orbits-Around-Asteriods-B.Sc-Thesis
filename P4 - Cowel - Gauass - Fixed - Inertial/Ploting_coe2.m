function Ploting_coe2(coe,T,color)

coe(:,3:6) = 180/pi*coe(:,3:6);
T = T/86400;

figure(2);
plot(T,coe(:,2),color);
xlabel('Days');
ylabel('Eccentricity');
hold on
legend('Cowel method','Slowly rotating asteroid method','Fixed-body method')

figure(3);
plot(T,mod(real(coe(:,3)),360),color);
xlabel('Days');
ylabel('Right ascension(degree)');
hold on
legend('Cowel method','Slowly rotating asteroid method','Fixed-body method')

figure(4);
plot(T,mod(real(coe(:,4)),360),color);
xlabel('Days');
ylabel('Inclination(degree)');
hold on
legend('Cowel method','Slowly rotating asteroid method','Fixed-body method')


end
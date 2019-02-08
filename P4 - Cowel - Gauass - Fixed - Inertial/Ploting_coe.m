function Ploting_coe(coe,T,color)

coe(:,3:6) = 180/pi*coe(:,3:6);
T = T/86400;

subplot(2,3,1);plot(T,coe(:,1),color);
xlabel('Days');
ylabel('Semimajor axis(km)');
hold on

subplot(2,3,2);plot(T,coe(:,2),color);
xlabel('Days');
ylabel('Eccentricity');
hold on

subplot(2,3,3);plot(T,mod(real(coe(:,3)),360),color);
xlabel('Days');
ylabel('Right ascension(degree)');
hold on

subplot(2,3,4);plot(T,mod(real(coe(:,4)),360),color);
xlabel('Days');
ylabel('Inclination(degree)');
hold on

subplot(2,3,5);plot(T,mod(real(coe(:,5)),360),color);
xlabel('Days');
ylabel('Argument of perigee(degree)');
hold on

subplot(2,3,6);plot(T,mod(real(coe(:,6)+coe(:,5)),360),color);
xlabel('Days');
ylabel('U (degree)');
hold on

end
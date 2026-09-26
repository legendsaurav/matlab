s=tf("s");
G=10*((s+2.5)^2)/((s+10)*(s^2 + 0.12));
figure;
bode(G);
grid on;









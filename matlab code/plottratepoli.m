%PLOT
figure(1);
hold on;
plot(t,I_vec, 'Color', "#D95319", "LineWidth",1) 
plot(t,B_vec, 'Color', "#EDB120", "LineWidth",1)
plot(t,D_vec, 'Color', "#7E2F8E", "LineWidth",1)
plot(t,R_vec, 'Color', "#77AC30", "LineWidth",1)
grid on;
xlabel('time')
ylabel('model variables')
title('EPG rate model Policarpine');
legend('I', 'B', 'D', 'R')

%PLOT R vs B
figure(2);
hold off;
plot(B_vec,R_vec, 'Color', "#D95319", "LineWidth",1) 
xlabel('B distruption')
ylabel('Remodelling')
title('B vs R');
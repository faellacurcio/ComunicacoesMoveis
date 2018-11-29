% clc 
% 
% powerMacro = 40-30; %(dbm)
% powerPico = 7-30; %(dbm)
% powerNoise = -115 - 30; % (db)
% T = 0.001; %(s)
% n = 4;
% d0 = 1000; %(m)
% alpha = 0.9;
% beta = 1.1; 
% B = 5*10^6; %(Hz)
% 
% TM_table = [0 1000; 1700 100];
% 
% potencia_recebida1 = powerMacro - pathLoss(hypot(TM_table(1,1),TM_table(1,2)));
% potencia_interferente1 = powerPico - pathLoss(hypot(TM_table(1,1)-1700,TM_table(1,2)-0));
% SINR1 = potencia_recebida1 - lin2db( db2lin(potencia_interferente1) + db2lin(powerNoise));
% 
% disp('capacity')
% capacity(1,SINR1)
% 
% % perda = db2lin(38.3 + 10*4*log10(300));
% % ruido = dbm2lin(-115);
% % SNR = db2lin(0);
% % lin2db(perda*ruido*SNR)
% 
% function C = capacity(uovertT, SINR) %(bits/s)
%     
% 	B = 5000000;
%     alpha = 0.9;
%     beta = 1.1;
%     C = min((uovertT)*alpha*B*log2(1+db2lin(SINR)/beta),10^7);
% 
% end
% 
% function pathloss = pathLoss(d)
% 
%     pathloss = 38.4 + 10*4*log10(d/1);
% 
% end
% 
% function result = normalizeData(data)
%     result = (data - min(data)) / ( max(data) - min(data) );
% end
% 
% function result = dbm2lin(input)
%     result = 10.^((input-30)/10);
% end
% 
% function result = db2lin(input)
%     result = 10.^((input)/10);
% end
% 
% function result = lin2dbm(input)
%     result = 10*log10(input) + 30;
% end
% 
% function result = lin2db(input)
%     result = 10*log10(input);
% end
TM_table = [4 3]

scatter(0,0,'*');
hold on; 

scatter(TM_table(1,1),TM_table(1,2),'o');

plot(...
    [0 cosd(atand((TM_table(1,2)/TM_table(1,1))))*hypot(TM_table(1,1),TM_table(1,2))], ...
    [0 sind(atand((TM_table(1,2)/TM_table(1,1))))*hypot(TM_table(1,1),TM_table(1,2))] ...
)

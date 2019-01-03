cd('C:\Users\Lior\master\evapocalc');
eilat = csvread('Eilat.csv', 1);
eilat = eilat(:,1:2);

Sedom = csvread('sedom.csv', 1);
Sedom = Sedom(:,1:2);

station = eilat;
stationDate = station(:,1);
stationRain = station(:,2);

intervals = stationDate(2:end)-stationDate(1:end-1);
stationDate(:,2) = [0;intervals];
stationDate(:,3:8) = datevec(stationDate(:,1));
intervalSeries = emprand(stationDate(:,2),10000,1);
RainAmountSeries = emprand(stationRain,10000,1);


[d_cdf,dx_cdf] = ecdf(stationRain);
csvwrite('DB\AmountCDF.csv', [d_cdf,dx_cdf]); 
[i_cdf,x_cdf] = ecdf(intervals);
csvwrite('DB\IntervalCDF.csv', [d_cdf,dx_cdf]); 

figure; plot(dx_cdf, d_cdf); hold on

figure; plot(x_cdf,i_cdf); hold on


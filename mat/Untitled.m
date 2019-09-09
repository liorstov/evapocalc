cd('C:\Users\Lior\master\evapocalc');
eilat = csvread('Eilat.csv', 1);
eilat = eilat(:,1:2);

Sedom = csvread('C:\Users\Lior\master\evapocalc\evapostats\sedom.csv', 1);
Sedom = Sedom(:,1:2);

station = Sedom;
stationDate = station(:,1);
stationRain = station(:,2);

SerialRain(:,1) = stationDate;
SerialRain(:,2) = stationRain;

rain = transpose([stationDate(1):stationDate(end)]);

%getting wet days in 6th column:
SerialRain(SerialRain(:,2)>0,3)=1;
%getting previous Wet/dry in 7
SerialRain(1,3)=0; %previous assumed dry
for i=2:length(SerialRain)
    SerialRain(i,3)=SerialRain(i,1) - SerialRain(i-1,1);
end

bla = SerialRain(SerialRain(:,3) == 1,:);

counter = 1;
for i=2:length(bla)
    if (bla(i,1)-bla(i-1,1) ~= 1)
        bla(i-1,4) = counter;
        counter = 1;
    else
        counter = counter + 1;
    end
    
end
bla = bla(bla(:,4) ~= 0,:);
eventDurationMean =  mean(bla(:,4));
dailyMean = mean(SerialRain(:,2));

% PWAW in 8 PWAD in 9
SerialRain(SerialRain(:,6)==1 & SerialRain(:,7)==1 ,8)=1;
SerialRain(SerialRain(:,6)==1 & SerialRain(:,7)==0 ,9)=1;
SerialRain(SerialRain(:,6)==0 ,8:9)=0; 
rain = transpose( [stationDate(1):stationDate(end)])

%the distribution is of the serial numbers of WAW\WAD:
DistWetDayPreviosWet=SerialRain(SerialRain(:,3)==1,2)
DistWetDayPreviosDry=SerialRain(SerialRain(:,3)>1,2)
DistWetDaySerial=SerialRain(:,2)

HstWetAfterWet = hist(DistWetDayPreviosWet,1:1:365);
HstWetAfterDry = hist(DistWetDayPreviosDry,1:1:365);
HstWet = hist(DistWetDaySerial,1:1:365);
PWetAfterWet(1) = 0;
PWetAfterDry(1) = 0;


DataLength= 59; %insert number of years in data - make sure data is made of whole seasons so no bias for months that are represented better than others. 
for i=2:365
    PWetAfterWet(i) = DistWetDayPreviosWet(i)/DistWetDaySerial(i-1);
    PWetAfterDry(i) = DistWetDayPreviosDry(i)/(DataLength-DistWetDaySerial(i-1)); 
end


%make sure no nans (check - shouldn't be any inf also):
PWetAfterWet(isnan(PWetAfterWet))=0
PWetAfterDry(isnan(PWetAfterDry))=0

P_WAW=PWetAfterWet; P_WAD=PWetAfterDry % that is the empiric distribution it is more resonable to smooth % it or to fit a gaussian:

%%%%%%%%%%%%%%%%% running window for wet day serial dist
%emp hist - running window of 'WindowSize' days

WindowSize = 15;
PWetAfterWetSmooth = conv(PWetAfterWet,ones(WindowSize,1)/WindowSize,'same');
PWetAfterDrySmooth = conv(PWetAfterDry,ones(WindowSize,1)/WindowSize,'same');

P_WAW=PWetAfterWetSmooth; P_WAD=PWetAfterDrySmooth;

[d_cdf,dx_cdf] = ecdf(stationRain);
csvwrite('DB\AmountCDF.csv', [d_cdf,dx_cdf]); 
[i_cdf,x_cdf] = ecdf(intervals);
csvwrite('DB\IntervalCDF.csv', [d_cdf,dx_cdf]); 

figure; plot(dx_cdf, d_cdf); hold on

figure; plot(x_cdf,i_cdf); hold on


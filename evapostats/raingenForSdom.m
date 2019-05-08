
%A=0.8502*(mean daily depth)-0.08192
A=0.8502*(2.3)-0.08192
B=0.76;
%Fd=0.06441(mean events in year)-0.2562
Fd=0.06441*(9)-0.2562;

%Fw=-1.897(mean events durations)(-1.187)+1.896
Fw=-1.897*(1.3)^(-1.187)+1.896;

%get rain
cd('C:\Users\Lior\master\evapocalc\evapostats');
eilat = csvread('Eilat.csv', 1);
eilat = eilat(:,1:2);

Sedom = csvread('C:\Users\Lior\master\evapocalc\evapostats\sedom.csv', 1);
Sedom = Sedom(:,1:2);

station = eilat;
stationDate = station(:,1);
stationRain = station(:,2);

pd = fitdist(stationRain,'wbl')
xx = 1:1:365
yy=wblpdf(xx,pd.A,pd.B); %fitted line
plot(yy)
%data = RainGen(DataLength,A,B,DepthLimit,Fw*P_WAW,Fd*P_WAD);

cd('C:\Users\Lior\master\evapocalc\evapostats');
eilat = csvread('Eilat.csv', 1);
eilat = eilat(:,1:2);

Sedom = csvread('C:\Users\Lior\master\evapocalc\evapostats\sedom.csv', 1);
Sedom = Sedom(:,1:2);

station = eilat;
stationDate = station(:,1);
stationRain = station(:,2);

SerialRain(:,1) = stationDate;
SerialRain(:,2) = stationRain;

rain = transpose([stationDate(1):stationDate(end)]);

%getting wet days in 6th column:
SerialRain(SerialRain(:,2)>0,3)=1;
%getting previous Wet/dry in 7
SerialRain(1,3)=0; %previous assumed dry

SerialRain(:,4) =  year(datetime(SerialRain(:,1), 'ConvertFrom','datenum'));
[a,b] = histc(SerialRain(:,4),unique(SerialRain(:,4)));
meanEventsPerYear = mean(a);
for i=2:length(SerialRain)
    SerialRain(i,3)=SerialRain(i,1) - SerialRain(i-1,1);
end

eventDurationMean =  mean(SerialRain(:,3));
dailyMean = mean(SerialRain(:,2));
meanAnualRain = mean(splitapply(@sum, SerialRain(:,2), findgroups(SerialRain(:,4))));

%A=0.8502*(mean daily depth)-0.08192
pd.A = 0.8502*(dailyMean)-0.08192;
pd.B = 0.76;
%Fd=0.06441(mean events in year)-0.2562
Fd=0.06441*(meanEventsPerYear)-0.2562;

%Fw=-1.897(mean events durations)(-1.187)+1.896
Fw=-1.897*(eventDurationMean)^(-1.187)+1.896;
%P_WAD  with gauss!
 a1=0.1748; b1=140.2; c1=79.33;
 x=1:1:365;
 P_WAD=a1*exp(-((x-b1)/c1).^2);
 %P_WAW  with gauss!
 a1=0.5828; b1=142.8; c1=96.93;
 P_WAW = a1*exp(-((x-b1)/c1).^2);
 
data = RainGen(1000,pd.A,pd.B,0.1,Fw*P_WAW,Fd*P_WAD);

csvwrite('..\DB\RainSeriesEilat.csv',data)   
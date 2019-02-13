function[SyntRainsData] = RainGen(SampleLength,A,B,DepthLimit,P_WAW,P_WAD)

%generate rains
SyntRainsData=[];
for y = 1:SampleLength
    %initial matrix for year
    %column for the year number
    SyntRain(1:365,1) = y;
    %column for the rainfall depth  - first get random numbers, next fill
                                        %with 1/0 according to the WAW/WAD probabilities
    SyntRain(1:365,3)=rand(365,1);
    
    %get random numbers for deriving rain depth from the rainfall
    %distribution
    seed = rand(365,1);
    
    %first day follow a dry day
    %column 2 is for the serial day of the year
    SyntRain(1,2) = 1;
    % high probabilities of wet days result with many wet day occurances,
    % low probabilities result with fewer:
    if SyntRain(1,3) <= P_WAD(1)
        %if WetDay, choose rain depth randomly from prob function 
        SyntRain(1,3) = wblinv(seed(1),A,B)+DepthLimit;
    else                                                    
        SyntRain(1,3)=0;%if NO>>DryDay 
    end
    %repeat for all days of the year
    for i=2:365
        SyntRain(i,2) = i;
        %choose a distribution conditiond on wet/dry previous day
        if  SyntRain(i-1,3)>0
            WetDayProb = P_WAW(i);
        else
            WetDayProb = P_WAD(i);
        end
  
        if SyntRain(i,3)<=WetDayProb
            SyntRain(i,3) = wblinv(seed(i),A,B)+DepthLimit ; 
        else                                                    
            SyntRain(i,3)=0;%if NO>>DryDay 
        end
    end
    SyntRainsData(length(SyntRainsData)+1:length(SyntRainsData)+length(SyntRain),:)=SyntRain;
end

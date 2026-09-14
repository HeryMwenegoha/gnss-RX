% By Hery A Mwenegoha (C) 2019 
% Read Rinex 3.02 File from CDDIS NASA Database
% Return a structure to each SV with ephemeris parameters GPS Only
% Reads the Broadcast Ephemeris File Unzipped CDDIS Files
% Files under IGS Folder, Otherwise direct to Database

% Assume DateTime of Flight:-
% dateTime_of_flight_yyyy_mm_dd_hh_min_sec   = [2019 03 02 14 00 00]
% dateTime_of_flight_yyyy_GPSDD_hh_min_sec   = [2019 61 14 00 00]
function [SV, Iono] = readEphemeris(filename)
fileID = fopen(filename);

% fgets reads new line including carriage return
% fgetl reads line excluding new line character

epoch = 0;
SvNum = 0;
IDs   = [];
flag_start = false;

SV=struct('ID',cell(1,32));
for iLine = 1:1000
    line_ex = fgetl(fileID);
    
    if line_ex == -1
        fclose(fileID);
        break;
    end
    
    if isempty(line_ex) || contains(line_ex, '****') 
        continue;
    end
    
    if contains(line_ex, 'GPSA') && iLine == 3
        cc     = strsplit(line_ex, ' '); 
        alpha0 = str2num(cc{2});
        alpha1 = str2num(cc{3});
        alpha2 = str2num(cc{4});
        alpha3 = str2num(cc{5}); 
        Iono.alpha = [alpha0, alpha1, alpha2, alpha3];
    elseif contains(line_ex, 'GPSB') && iLine == 4
        cc     = strsplit(line_ex, ' ');
        beta0  = str2num(cc{2});
        beta1  = str2num(cc{3});
        beta2  = str2num(cc{4});
        beta3  = str2num(cc{5}); 
        Iono.beta = [beta0, beta1, beta2, beta3];
    end
    
    if contains(line_ex, 'TIME SYSTEM CORR')
        TYPE    = line_ex(1:4);    % TYPE - [GAUT GPUT SBUT GLUT GPGA GLGP QZGP QZUT  BDUT]
        SKIP    = line_ex(5);
        a0      = line_ex(6:22);   % [sec]
        a1      = line_ex(23:38);  % [sec/sec]
        Tref_UTC= line_ex(39:45);  % [sec]
        Wref_UTC= line_ex(46:50);  % [GPS WEEK] 
        SKIP    = line_ex(51);
        SNN     = line_ex(52:56);
        SKIP    = line_ex(57);
        UTC     = line_ex(58:59);
        %[TYPE  SKIP a0 a1 Tref_UTC Wref_UTC SKIP SNN SKIP UTC]
    end
        
    % found end of header line, we expect GPS satellites
    if contains(line_ex,'END OF HEADER')
        flag_start = true;
        line_ex;
        continue;
    end
    
    if flag_start == false
        continue; % skip code until flag is changed indicating gps eph.
    end
    
    % from here we expect GPS Ephemeris messages : Data Stuff
    if contains(line_ex, 'G')
        % Update Space Vehicle Number
        SvNum = SvNum + 1;
        
        CONSTELLATION = line_ex(1);
        PRN           = str2double(line_ex(2:3));
        prn_          = sprintf('%s%u',CONSTELLATION, PRN);       
        
        % Time of Clock in GPS time
        TOC.yyyy = str2num(line_ex(4:8));
        TOC.mm   = str2num(line_ex(9:11));
        TOC.dd   = str2num(line_ex(12:14)); % Day on month
        TOC.hh   = str2num(line_ex(15:17));
        TOC.min  = str2num(line_ex(18:20));
        TOC.sec  = str2num(line_ex(21:23));
        
        toc_g_time=g_time(datetime(TOC.yyyy,...
                                   TOC.mm,...
                                   TOC.dd,...
                                   TOC.hh,...
                                   TOC.min,...
                                   TOC.sec));
                
        % TotalDays Since Beg.Of.Year Until Reference Epoch Day
        % Note : reference epoch has not finished
        % actualdays =  totaldays-1
        if TOC.mm > 1
            days_of_month = 28;      % account for 28 days of each month
            DayOfYear    = days_of_month * (TOC.mm -1);
            
            % days to adjust each month to account for 28daysPerMonth Assumption 
            add_days = [3 0 3 2 3 2 3 3 2 3 2 3];
        
            % Total number of Months to Adjust Days excluding the current month
            months_adj = (TOC.mm -1);
        
            % TotalDaysElapsed_In_This_Year = Current_Day_Of_This_Month + Adjusted_Days_of_Prev_Months
            DayOfYear = TOC.dd + days_of_month.*months_adj + sum(add_days(1:months_adj), 2);
        else
            DayOfYear = TOC.dd;  % If first month just get number of days passed. Jan 1st being day 1       
        end
        
        % Satellite Clock parameters
        sat.ao = str2num(line_ex(24:42));
        sat.a1 = str2num(line_ex(43:61)); %
        sat.a2 = str2num(line_ex(62:80));
        IDs(end+1) = PRN;
        if numel(find(IDs == PRN)) > 1
            fprintf('...duplicate PRN %u found ... \n',PRN);
        end
        
        % Start Packing Messages in Satellite Structure object
        % Reference epoch
        currentYYYY     = TOC.yyyy;
        YYYYref         = 1980;
        YearsFrom1980   = (currentYYYY-YYYYref);
        dayInThoseYears = YearsFrom1980*365;
        leapDays        = 0;
        for idx=YYYYref:currentYYYY
            if mod(idx,4) == 0 
                if mod(idx,100) == 0
                    %skip
                    if mod(idx,400) == 0
                        % put back in
                        leapDays = leapDays + 1; 
                    end
                else
                    leapDays = leapDays + 1; 
                end
            end
        end
        
        % Total days since Jan 6th 1980
        Jan5daysOf1980     = 5;  % remove the first 5days of January1980
        FulldaysElapsed    = (dayInThoseYears + leapDays + (DayOfYear-1)) - Jan5daysOf1980;
        FullSecondsElapsed = FulldaysElapsed*86400; % until last full day
        TotalWeeks         = floor(FullSecondsElapsed/604800);
        DaysIntoWeek       = floor((FullSecondsElapsed/604800 - TotalWeeks)*7);
        SecondsIntoWeek    = mod(FullSecondsElapsed,604800) + TOC.hh.*3600 + TOC.min.*60 + TOC.sec;
        
        % - ALERT TOC DOW - GPS
        prn         = PRN;
        SV(prn).ID  = prn;
        SV(prn).TOC = SecondsIntoWeek;  
        SV(prn).TOC_g_time = toc_g_time;
        SV(prn).A0  = sat.ao;
        SV(prn).A1  = sat.a1;
        SV(prn).A2  = sat.a2;        
        continue;
    end
    
    % 7 lines then reset epoch
    epoch = epoch + 1;
    switch(epoch)
        case 1
            EMPTY  = line_ex(1:4);
            IODE   = line_ex(5:23); 
            CRS    = line_ex(24:42); % metres
            DELTA_N= line_ex(43:61); % rad/sec
            M0     = line_ex(62:80); % rad
            
            SV(prn).IODE     = str2num(IODE);
            SV(prn).CRS      = str2num(CRS);
            SV(prn).DELTA_N  = str2num(DELTA_N);
            SV(prn).M0       = str2num(M0);
            
        case 2
            EMPTY  = line_ex(1:4);
            CUC    = line_ex(5:23);  % rad - ArgofLat Cosine Harmonic Corr
            e0     = line_ex(24:42); % Eccentricity
            CUS    = line_ex(43:61); % rad
            SQRT_A = line_ex(62:80); % sqrt(m)
            
            SV(prn).CUC     = str2num(CUC); 
            SV(prn).e0      = str2num(e0);
            SV(prn).CUS     = str2num(CUS);
            SV(prn).SQRT_A  = str2num(SQRT_A);

        case 3
            EMPTY  = line_ex(1:4);
            TOE    = line_ex(5:23);  % sec  - of GPS week
            CIC    = line_ex(24:42); % rad
            OMEGA0 = line_ex(43:61); % rad
            CIS    = line_ex(62:80); % rad 
            
            SV(prn).TOE     = str2num(TOE); 
            SV(prn).CIC     = str2num(CIC);
            SV(prn).OMEGA0  = str2num(OMEGA0);
            SV(prn).CIS     = str2num(CIS);

        case 4
            EMPTY   = line_ex(1:4);
            IO      = line_ex(5:23);  % rad  - inclination
            CRC     = line_ex(24:42); % metres
            omega   = line_ex(43:61); % rad
            OMEGADOT= line_ex(62:80); % rad/sec 
            
            SV(prn).IO        = str2num(IO); 
            SV(prn).CRC       = str2num(CRC);
            SV(prn).omega     = str2num(omega);
            SV(prn).OMEGADOT  = str2num(OMEGADOT);

        case 5
            EMPTY   = line_ex(1:4);
            IDOT    = line_ex(5:23);    % rad/s  - inclination
            CODESL2 = line_ex(24:42);   % []
            GPSWEEK = line_ex(43:61);   % Continous GPS week number not mod(1024) [] - correct
            L2PDATAFLAG=line_ex(62:80); % [] 
            
            SV(prn).IDOT = str2num(IDOT); 
            SV(prn).CODESL2  = str2num(CODESL2);
            SV(prn).GPSWEEK  = str2num(GPSWEEK);
            SV(prn).L2PDATAFLAG  = str2num(L2PDATAFLAG);
            
        case 6
            EMPTY     = line_ex(1:4);
            SVACCURACY= line_ex(5:23);  % [metres]
            SVHEALTH  = line_ex(24:42); % Bits 17-22
            TGD       = line_ex(43:61); % [seconds] Group delay 
            IODC      = line_ex(62:80); % Issue of Data Clock []  
            
            SV(prn).ACCURACY = str2num(SVACCURACY); % SV   Signal In Space Accuracy [SISA] ACCURACY [m]
            SV(prn).HEALTH   = str2num(SVHEALTH);   % SV   HEALTH
            SV(prn).TGD      = str2num(TGD);        % TGD  
            SV(prn).IODC     = str2num(IODC);       % IODC 
            %str2num(TGD)*3e8
            
        case 7
            EMPTY     = line_ex(1:4);
            ZCount    = line_ex(5:23);  % [sec] of GPS Week - Also transmission time of message in Zcount
            FIT       = line_ex(24:42); % [hours] ICD-GPS-200,20.3.4.4; 0=4hrs, 1=6hrs
            SV(prn).ZCount = str2num(ZCount); 
            SV(prn).FIT    = str2num(FIT);
            epoch = 0;
    end
end

% If fileId is open then close it
if any(fopen('all') == fileID)
    fclose(fileID);
end
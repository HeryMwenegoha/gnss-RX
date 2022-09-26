classdef g_time 
    % Can be used to hold any time object in calendar format
    properties
        date_time;
        fractional_second;
    end
    
    methods
        function obj = g_time(datetime, fractional_second)
            if ~isa(datetime, 'datetime')
                error('g_time::g_time input not datetime');
            end
            obj.date_time         = datetime;
            if nargin > 1
                obj.fractional_second = fractional_second;
            else
                obj.fractional_second = 0;
            end
        end
        
        %{
        function obj3 = plus(obj1, obj2)
            sum_fractional =   obj1.fractional_second +...
                               obj2.fractional_second;
            dt = obj1.date_time +...
                 obj2.date_time +...
                 seconds(floor(sum_fractional));
             ff = sum_fractional - floor(sum_fractional);
            obj3 = g_time(dt, ff);
        end
        %}
        
        function obj3 = plus(obj1, obj2_duration)
            if ~isa(obj2_duration, 'duration')
                error('g_time::plus::Please specify duration object');
            end
             obj2_duration.Format = 's';
             obj2_dur_ss  = obj2_duration;
             new_datetime = obj1.date_time +  obj2_dur_ss;
             new_ff       = obj1.fractional_second;
             obj3         = g_time(new_datetime, new_ff);
        end
        
        
        function y = colon(a, d, b)
            obj1 = a;
            aa   = obj1.date_time + seconds(obj1.fractional_second);
            if nargin == 2
                obj2 = d;
            elseif nargin == 3
                obj2 = b;
            end
            bb = obj2.date_time + seconds(obj2.fractional_second);
            
            if nargin == 2
                y  = aa:bb;
            elseif nargin == 3
                if isa(d, 'duration') == true
                    y  = aa:d:bb;
                elseif isa(d, 'double') == true
                    y  = aa:seconds(d):bb;
                end
            end
        end
        
        
        function y = minus(obj1, obj2)
            diff_fractional = obj1.fractional_second -...
                              obj2.fractional_second;
            dt        = obj1.date_time - obj2.date_time;% duration object
            ff        = diff_fractional;                % could be negative
            z         = dt + seconds(ff);
            y.seconds = seconds(z);
            y.days    = days(z);
            y.hour    = hours(z);
            y.minute  = minutes(z);
        end
        
        function y=g_week(obj)
            y=g_time.gps_week(obj);
        end
        
        function y=g_sow(obj)
            y=g_time.gps_sow(obj);
        end
        
        function y=g_doy(obj)
            y=g_time.gps_doy(obj);
        end
        
        function y=g_diw(obj)
            y=g_time.gps_diw(obj);
        end        
    end
    
    methods(Static)
        % can be accessed statically
        function y = formatter(obj)
           if isa(obj, 'g_time') == true
                % do nothing
            elseif isa(obj, 'datetime') == true 
                % format convertion
                obj = g_time(obj, 0);
            else
                warning("g_time::formatter input not handled");
                y = [];
                return;
           end 
            y = obj;
        end
        
        function y = gps_week(obj)
            gps_g_time  = g_time(datetime(1980,01,06,00,00,00), 0);
            obj         = g_time.formatter(obj);
            diff_g_time = obj - gps_g_time;
            y=floor(diff_g_time.seconds/604800);
        end
        
        
        function y = gps_sow(obj)
           obj = g_time.formatter(obj);
           y   =  g_time.gps_diw(obj)*86400 +...
                  obj.date_time.Hour*3600   + ...
                  obj.date_time.Minute*60   + ...
                  obj.date_time.Second      + ...
                  obj.fractional_second;          
        end
        
        
        function y = gps_doy(obj)
             obj = g_time.formatter(obj);
             yyyy        = obj.date_time.Year;
             diff_g_time = obj - g_time(datetime(yyyy,01,01,00,00,00),0);
             y           = floor(diff_g_time.days) + 1;
        end
        
        
        function y = gps_diw(obj)
            % Carefull with the interpretation of DiW-Days into Week
            % GPS Has 7 days - 1 full day is equivalent to the passage of
            % 24 hours or 86400seconds
            % day 1     - Sunday 
            % day 2     - Monday
            % day 3     - Tuesday
            % day 4     - Wednesday
            %      .5   - Thursday-Midday
            % day 5     - Thursday
            % day 6     - Friday
            %      .5   - Saturday-Midday
            % day 7     - Saturday
            % GPS days into week i.e. number of days that have fully passed
            % If diw = 1, means 1 full day has passed 
            % If diw = 4, means 4 full days have passed
            gps_g_time = g_time(datetime(1980,01,06,00,00,00), 0);
            obj = g_time.formatter(obj);
            diff_g_time = obj - gps_g_time;
            y=mod(floor(diff_g_time.days),7);
        end
        
    end
end


% The current day of the year
%{
 if isa(obj, 'g_time') == true
     % Get the year of the passed object
     yyyy        = obj.date_time.Year;
     diff_g_time = obj - g_time(datetime(yyyy,01,01,00,00,00),0);
     y           = floor(diff_g_time.days) + 1;
elseif isa(obj, 'datetime') == true
     yyyy        = obj.Year;
     diff_g_time = obj - datetime(yyyy,01,01,00,00,00);
     y           = floor(diff_g_time.days) + 1;              
else
    warning("g_time::gps_week input not handled");
 end 
%}
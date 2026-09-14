% Hery A Mwenegoha copyright (c) 2020 - 2026

classdef gTime 
    properties
        t = [];
        f = [];
    end
    
    methods
        function obj=gTime(t,f_)
            % Declare a gTime object to handle separate integer and
            % fractional time objects. 
            f=0;
            if nargin == 1
            else
                f = f_; % Fractional part
            end
            obj.t = t + fix(f);
            obj.f = f - fix(f);
        end
        
        function r = plus(obj1, obj2)
            tt = (obj1.t + obj2.t);
            ff = (obj1.f + obj2.f);            
            r  = gTime(tt,ff);
        end
        
        function r = minus(obj1, obj2)
            tt = (obj1.t - obj2.t);
            ff = (obj1.f - obj2.f);
            r  = gTime(tt,ff);
        end
        
        function r = le(obj,b)
            if isa(b, 'gTime')
                b = b.t + b.f;
            end
            
            if isa(obj, 'gTime')
                if (obj.t + obj.f) <= b
                    r = true;
                else
                    r = false;
                end
            else
                error('gTime:lt obj not gTime object');
            end
        end
 
        function r = lt(obj,b)
            if isa(b, 'gTime')
                b = b.t + b.f;
            end
            
            if isa(obj, 'gTime')
                if (obj.t + obj.f) < b
                    r = true;
                else
                    r = false;
                end
            else
                error('gTime:lt obj not gTime object');
            end
        end
        
        function r = ge(obj,b)
            if isa(b, 'gTime')
                b = b.t + b.f;
            end
            
            if isa(obj, 'gTime')
                if (obj.t + obj.f) >= b
                    r = true;
                else
                    r = false;
                end
            else
                error('gTime:lt obj not gTime object');
            end
        end
        
        function r = gt(obj,b)
            if isa(b, 'gTime')
                b = b.t + b.f;
            end
            
            if isa(obj, 'gTime')
                if (obj.t + obj.f) > b
                    r = true;
                else
                    r = false;
                end
            else
                error('gTime:lt obj not gTime object');
            end
        end
        
        
        function r = times(obj,b)
            if isa(obj, 'gTime') && isa(b, 'double')
                r=obj.t*b + obj.f*b;
            elseif isa(obj, 'double') && isa(b, 'gTime')
                r=obj*b.t + obj*b.f;
            else
                error('gTime::times::not handled');
            end
        end
        
        
        function r = power(obj,b)
            if isa(obj, 'gTime') && isa(b, 'double')
                r = (obj.t + obj.f).^b;
            elseif isa(obj, 'double') && isa(b, 'gTime')
                error('gTime::power::not handled');
            else
                error('gTime::power::not handled');
            end
        end

        function r = mpower(obj,b)
            if isa(obj, 'gTime') && isa(b, 'double')
                r = (obj.t + obj.f)^b;
            elseif isa(obj, 'double') && isa(b, 'gTime')
                error('gTime::mpower::not handled');
            else
                error('gTime::mpower::not handled');
            end
        end
        
        function r = mtimes(obj,b)
            if isa(obj, 'gTime') && isa(b, 'double')
                r=obj.t*b + obj.f*b;
            elseif isa(obj, 'double') && isa(b, 'gTime')
                r=obj*b.t + obj*b.f;
            else
                error('gTime::mtimes::not handled');
            end
        end
    end
    
    methods(Static)
        function [I,F] = formatter(input)
            I=round(input);
            F=input-I;
        end
    end
end
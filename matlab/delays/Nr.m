% By Hery A Mwenegoha (C) 2019

classdef Nr < handle
    properties(Access = private)
        N_i;
        N_i_old;
    end
    
    properties(Constant)
       P=[27.9947636532892,   -2.393938120604054
          28.50742139015308,  -2.9077950455109636
          28.994453351261892, -3.40151466177079
          29.50712401737691,  -3.9254725641290853
          30.000594745552334, -4.449478951179147
          30.506917149358706, -5.013856924906909
          30.981187939594538, -5.5379117966487375
          31.51316197765825,  -6.1426290339263545
          31.987406909391808, -6.646481950765412];
 
        p1  = -1.071  %(-1.092, -1.05)
        p2  = 27.62   %(27, 28.25)
    end
    
    methods
        function this = Nr
            this.N_i     = randi([-20,20],1);
            this.N_i_old = this.N_i;
        end 
        
 
        function y = N(this)
            % Case slip first
            if randi([1,1000]) == 400
                this.N_i = randi([-50,50],1);
            end
            y = this.N_i;
        end
        
        function y =detect(this)
            % Try to detetc
            if this.N_i_old == this.N_i
                y = false;
            else
                this.N_i_old = this.N_i;
                y          = true;
            end
        end
    end
    
    methods(Static)
        function  y = slip(El_in_deg)
            % get a random number that is normally distributed
            v = randn(1000,1);
            
            % CNo based on elevation
            cno = CNO.cno(El_in_deg);
            
            % probability of cycle-slip
            f   = Nr.prob_log10(cno);
            
            % actual probability of getting a slip
            p    = 10.^f;
            
            % probability of not getting a slip
            pinv  = (1-p) + p/2;
            
            iip   = v >  norminv(pinv);
            iip2  = v < -norminv(pinv);
            
            % z-score of pinv
            if any(iip)
                slip_happened = true;
            elseif any(iip2)
                slip_happened = true;
            else
                slip_happened = false;
            end
            
            y = slip_happened;
        end
        
        % probability-log10 of a cycle slip
        function f=prob_log10(cno)
            ip=cno > 32;
            cno(ip)=32;
            
            ip=cno < 26.1;
            cno(ip)=26.1;
            
            f = Nr.p1*cno + Nr.p2;
        end
        
        function plot(cno)
            plot(Nr.P(:,1), Nr.P(:,2));
            plot(cno, Nr.prob_log10(cno));
        end
    end
    
end
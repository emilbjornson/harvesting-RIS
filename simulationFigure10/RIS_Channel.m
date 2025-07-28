classdef RIS_Channel< handle
    properties
        c = 299792458; % m/s
        freq = 2.4e9; % Hz
        lambda = [];
        %%
        Polarization = true;
        %% TX
        TX_pos = [0;-3;4];
        G_TX_dBi = 0;%15.478906952036844;
        G_TX = [];
        TX_polarization = [1;0;0];
        P_TX = 1; % Watt
        %% RIS
        D_RIS = [];
        d_RIS =  0.5; % times lambda
        N_RIS_row = 10;
        N_RIS_col = 10;
        RIS_elements_pos = {}; % It is located in XY plane
        RIS_center_pos = [0;0;0];
        RIS_polarization = [1;0;0]; % It is in X direction
        Phase_shift_set = [0;pi/2;pi;3*pi/2];
        Optimal_Phase_Shifts = [];
        Sub_Optimal_Phase_Shifts = [];
        %% RIS Channel
        G_RIS_For_TX = [];
        G_RIS_For_EHU = [];
        Time_delay_TX_RIS = [];
        Time_delay_RIS_EHU = [];
        phi_h = [];
        phi_g = [];
        h = [];
        g = [];
        z = [];
        %% Energy Harvesting Unit (EHU)
        EHU_pos = [0;1;2];
        G_EHU_dBi = 0;%3.267771379877908;
        G_EHU = [];
        EHU_polarization = [1;0;0];
        %% Harvested power
        Total_incident_Power = [];
        Total_received_Power = [];
        Total_harvested_Power = [];
    end
    methods
        function obj = RIS_Channel(varargin)
            if nargin>0
                obj.N_RIS_row = varargin{1};
                obj.N_RIS_col = varargin{2};
            end
            obj.lambda = obj.c/obj.freq;
            obj.G_RIS_For_TX = zeros(obj.N_RIS_row, obj.N_RIS_col); % Should be modified.
            obj.RIS_elements_pos = cell(obj.N_RIS_row, obj.N_RIS_col);
            obj.Time_delay_TX_RIS = zeros(obj.N_RIS_row, obj.N_RIS_col);
            obj.Time_delay_RIS_EHU = zeros(obj.N_RIS_row, obj.N_RIS_col);
            obj.D_RIS = obj.lambda*obj.N_RIS_row*obj.d_RIS;
            obj.G_TX = 10^(obj.G_TX_dBi/10);
            obj.G_EHU = 10^(obj.G_EHU_dBi/10);
        end
        function  obj = Harvested_power_func(obj)
            Incident_power = zeros(obj.N_RIS_row,obj.N_RIS_col);
            RX_power = zeros(obj.N_RIS_row,obj.N_RIS_col);
            for i = 1:obj.N_RIS_row
                for j = 1:obj.N_RIS_col
                    obj.RIS_elements_pos{i,j} = obj.lambda*[(i-(obj.N_RIS_row+1)/2)*obj.d_RIS;(j-(obj.N_RIS_col+1)/2)*obj.d_RIS;0];
                    obj.Time_delay_TX_RIS(i,j)=norm(obj.TX_pos-obj.RIS_elements_pos{i,j})/obj.c;
                    obj.Time_delay_RIS_EHU(i,j)=norm(obj.EHU_pos-obj.RIS_elements_pos{i,j})/obj.c;
                    if obj.Polarization
                        roh_n_t = obj.TX_polarization-dot((obj.RIS_elements_pos{i,j}-obj.TX_pos),obj.TX_polarization)*(obj.RIS_elements_pos{i,j}-obj.TX_pos)/(norm((obj.RIS_elements_pos{i,j}-obj.TX_pos))).^2;
                        roh_n_t=roh_n_t/norm(roh_n_t);
                        roh_n_r = obj.EHU_polarization-dot((obj.RIS_elements_pos{i,j}-obj.EHU_pos),obj.EHU_polarization)*(obj.RIS_elements_pos{i,j}-obj.EHU_pos)/(norm((obj.RIS_elements_pos{i,j}-obj.EHU_pos))).^2;
                        roh_n_r=roh_n_r/norm(roh_n_r);
                    else
                        roh_n_t = obj.TX_polarization;
                        roh_n_r = obj.EHU_polarization;
                    end
                    A = obj.RIS_elements_pos{i,j}+obj.lambda*[obj.d_RIS/2;obj.d_RIS/2;0];
                    B = obj.RIS_elements_pos{i,j}+obj.lambda*[-obj.d_RIS/2;obj.d_RIS/2;0];
                    C = obj.RIS_elements_pos{i,j}+obj.lambda*[-obj.d_RIS/2;-obj.d_RIS/2;0];
                    D = obj.RIS_elements_pos{i,j}+obj.lambda*[obj.d_RIS/2;-obj.d_RIS/2;0];
                    %%  RIS Gain for TX
                    obj.G_RIS_For_TX(i,j) = (Spherical_triangle_Area(obj,A,B,C,obj.TX_pos)+Spherical_triangle_Area(obj,A,D,C,obj.TX_pos))*abs(dot(obj.RIS_polarization,roh_n_t))^2/4/pi;
                    %%  RIS Gain for EHU
                    obj.G_RIS_For_EHU(i,j) = (Spherical_triangle_Area(obj,A,B,C,obj.EHU_pos)+Spherical_triangle_Area(obj,A,D,C,obj.EHU_pos))*abs(dot(obj.RIS_polarization,roh_n_r))^2/4/pi;
                    %% Incident Power
                    Incident_power(i,j) = obj.P_TX*obj.G_TX*obj.G_RIS_For_TX(i,j);
                    %% Receive Power
                    RX_power(i,j) = Incident_power(i,j)*obj.G_EHU*obj.G_RIS_For_EHU(i,j);
                end
            end
            obj.Total_incident_Power = sum(Incident_power(:));
            obj.Total_received_Power = (sum(sqrt(RX_power(:)))).^2;
        end
        function  y = Spherical_triangle_Area(~,A,B,C,O)
            %% Inputs
            % O is the position vector of the center of the sphere
            % A, B, and C are the position vectors of the vertices of the triangle
            %% Output
            % y is the normalized area of the projection of the triangle on
            % the sphere surface
            OA = O-A;
            OB = O-B;
            OC = O-C;
            a =  abs(atan2(norm(cross(OB,OC)), dot(OB,OC)));
            b =  abs(atan2(norm(cross(OA,OC)), dot(OA,OC)));
            cc =  abs(atan2(norm(cross(OA,OB)), dot(OA,OB)));
            s = (a+b+cc)/2;
            y = 4*atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-cc)/2)));
        end
        function obj = RIS_CH(obj)
            obj.Harvested_power_func();
            tmp = obj.Time_delay_TX_RIS*obj.c/obj.lambda;
            obj.phi_h = (tmp - floor(tmp))*2*pi;
            tmp = obj.Time_delay_RIS_EHU*obj.c/obj.lambda;
            obj.phi_g = (tmp - floor(tmp))*2*pi;
            obj.h = sqrt(obj.G_RIS_For_TX).*exp(1j*obj.phi_h);
            obj.g = sqrt(obj.G_RIS_For_EHU).*exp(1j*obj.phi_g);
            obj.z = obj.h(:).*obj.g(:);
        end
        function y = Harvesting_efficiency(~,x,a,b,Psat)
            y = Psat*(1./(1+exp(-1*a*(x-b)))-1/(1+exp(a*b)))/(1-1/(1+exp(a*b)));
            y = y(:);
        end
        function [rx_power, inc_power] = Infinit_surface(obj,d)
            a = obj.TX_pos;
            b = obj.EHU_pos;
            if obj.Polarization
                fun_rx = @(x,y) sqrt(((y-a(2)).^2+a(3).^2).*(((y-b(2)).^2+b(3).^2)))./(((x-a(1)).^2 + (y-a(2)).^2+a(3).^2).*((x-b(1)).^2 + (y-b(2)).^2+b(3).^2)).^(5/4);
                fun_inc = @(x,y) ((((y-a(2)).^2+a(3).^2).^2)./((((y-a(2)).^2+a(3).^2).^2)+(x-a(1)).^2.*((y-a(2)).^2+a(3).^2)))./(((x-a(1)).^2 + (y-a(2)).^2+a(3).^2)).^(3/2);
            else
                fun_rx = @(x,y) 1./(((x-a(1)).^2 + (y-a(2)).^2+a(3).^2).*((x-b(1)).^2 + (y-b(2)).^2+b(3).^2)).^(3/4);
                fun_inc = @(x,y) 1./(((x-a(1)).^2 + (y-a(2)).^2+a(3).^2)).^(3/2);
            end
            tmp_inc = obj.P_TX.*a(3)/4/pi;
            inc_power = tmp_inc*(integral2(fun_inc,-d,d,-d,d));
            tmp_rx = obj.P_TX.*a(3)*b(3)/16/pi^2;
            rx_power = tmp_rx*(integral2(fun_rx,-d,d,-d,d))^2;
        end
    end
end
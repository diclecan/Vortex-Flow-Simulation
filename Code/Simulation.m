classdef Simulation
    %Wichtige funktionen und die Lösung der DGLs passiert hier

    properties %Die variablen die unter den funktionen weitergegeben werden. Diese müssen dann mit obj. erweitert werden
        nu, T,Nx, Ny, hx,hy, ht, Dx1, Dx2, Dy1, Dy2, L, R, P % wenn man obj in eine methode schreibt kann automaitsh per obj. auf die variable zugegriffen werden.
    end

    methods % Die wichtigen funktionen und die DGLs

        function obj = Simulation(N, nu, CFLx, u_max) % Start der Simulation bzw. grundeinstellungen    !!!u_max weg zum testen!!!
            obj.Nx = N; 
            obj.Ny = N; 
            obj.nu = nu; %grundlegendes festlegen

            obj.hx = (2*pi)/N; 
            obj.hy = (2*pi)/N; % Gitterweite
            obj.ht = (CFLx*obj.hx)/u_max;

            [obj.Dx1, obj.Dx2, obj.Dy1, obj.Dy2] = obj.createPeriodicDerivatives();%periodische ableitungsMatrizen erstellen
            A = obj.Dx2 + obj.Dy2; % Die poissonmatrix A*Psi = - omega, A wäre laplace
            [obj.L, obj.R, obj.P] = lu(A); %LR zerlegung von Poissonmatrix zum einfachen lösen. mit P*A = L*R --> P\(L*R) =A
        end
        % ab hier funktionen für die zeitschleife
        function Psi = PoissonSolver(obj,omega_vektor) %Poisson lösen
            
            zwischenvektor = (obj.L)\(-obj.P*omega_vektor);
            Psi = obj.R\zwischenvektor; % Psi kommt übrigens als vektor raus
        end

        function [u,v] = velocities(obj, Psi)%cauchy-Riemann lösen bzw einfach ableiten
            u = obj.Dy1*Psi; 
            v = -obj.Dx1*Psi;
        end

        function omega_vektor_next = wirbetransport(obj,u,v, omega_vektor)%Wirbeltransport glg lösen für Wirbelstärkeverteilung des nächsten Zeitschritts
            N_gesamt = obj.Nx * obj.Ny;
            u_matrix = spdiags(u, 0, N_gesamt, N_gesamt); 
            v_matrix = spdiags(v, 0, N_gesamt, N_gesamt); %MAcht geschwindigkeiten zu matruzen für linke seite des kommenden LGS
            diffusion = obj.nu*(obj.Dx2+obj.Dy2)*omega_vektor; %Explizite Diffusion
            
            %Schiefsymmetrie
            Conv_schiefsym = 0.5* ( (u_matrix*obj.Dx1+v_matrix*obj.Dy1) + (obj.Dx1*u_matrix+obj.Dy1*v_matrix) ); % 0.5 * (nicht konservative Konvektion + konservative Konvektion)

            A = speye(N_gesamt) + obj.ht * Conv_schiefsym; %RESULTIERENDE MATRIX
            omega_vektor_next = A\(omega_vektor + obj.ht*diffusion);

        end
        
        function [Dx1, Dx2, Dy1, Dy2] = createPeriodicDerivatives(obj)
            e = ones(obj.Nx, 1);
            Cx1 = spdiags([-e e], [-1 1], obj.Nx, obj.Nx);
            Cx1(1,end) = -1; Cx1(end,1) = 1; Cx1 = Cx1 / (2*obj.hx);
            Cx2 = spdiags([e -2*e e], [-1 0 1], obj.Nx, obj.Nx);
            Cx2(1,end) = 1; Cx2(end,1) = 1; Cx2 = Cx2 / (obj.hx^2);
            Ix = speye(obj.Nx);
            Dx1 = kron(Ix, Cx1); Dx2 = kron(Ix, Cx2);
            Dy1 = kron(Cx1, Ix); Dy2 = kron(Cx2, Ix);
        end
    end % ende der methoden


end % ende der klasse
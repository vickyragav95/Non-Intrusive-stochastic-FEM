clear;
clc;
% case number
q = 2;
% problem parameters
P = 3; % maximun gPC order
d = 1; % number of stochastic variable
M = factorial(P+d)/factorial(P)/factorial(d); %total number of terms in Pascal's simplex
kappa = 1E-2;

% beta(y) stochastic = 1+y^2
beta_cap = zeros(M,1);
beta_cap(1) = 1.5;
beta_cap(3) = 0.5;

% BCs
g1 = 1; 
g2 = 1;

nel = 50;
y_points = 50;

% num. of points
npoints = nel+1;

% domain (starting at x_A and of length L)
x_A = 0.0;
L = 1.0;
x_B = x_A+L;

% mesh spacing
delx = L/nel;

% mesh points
for ipoint = 1:npoints
  xmesh(ipoint) = x_A + (ipoint-1)*delx;
end

% mesh connectivity
for iel = 1:nel
  ien(iel,1) = iel;
  ien(iel,2) = iel+1;
end

% num. of (Gauss) quadrature points
nq = 2; % accurate up to degree 3 (e.g., linear basis and quadratic source term in: \int N_a s dx)
% local location (xi-coordinates) : quadrature points (interval [-1,1])
xiq = [-1.0/sqrt(3),1.0/sqrt(3)];
% weight : quadrature points
wq = [1.0,1.0];

% local shape function at quadrature points
% N(1) is left local shape function
% N(2) is right local shape function
% dN(1) is local (or xi) derivative of N1 
% dN(2) is local (or xi) derivative of N2
for iq = 1:nq 

  N(1,iq) = 0.5*(1-xiq(iq));
  N(2,iq) = 0.5*(1+xiq(iq));

  dN(1,iq) = -0.5;
  dN(2,iq) = 0.5;
  
end

[indx,Ckij] = lexico_table(P,d,1);
const = zeros(1,M);
const(1) = 1;% gPC coefficinets for const=1 [1 0 0 ... 0]

% [K_{AB}^{IJ}]{U_B^J}={F_A^I}
K_cell =  cell(M,M);
F_cell = cell(M,1);

for method=1:2
    for I=1:M
        for J=1:M
            % LHS matrix (dense storage)
            A = zeros(npoints,npoints); % eye(npoints,npoints);
            % RHS vector
            b = zeros(npoints,1);

            for iel = 1:nel
            %==============

              nshl = 2; % line or 1D element with linear basis
              Ae = zeros(nshl,nshl); % element-level LHS matrix
              be = zeros(nshl,1); %  element-level RHS vector    

              xl = [xmesh(ien(iel,1)),xmesh(ien(iel,2))]; % element-level nodes

              for iq = 1:nq
              %============

                %%%% prep data - begin
                Nq = N(:,iq); % local shape functions at this quadrature point
                dNq = dN(:,iq); % local (or xi) derivative of shape functions at this quadrature point

                % xglobal = 0.0; % global x location at this quadrature point - needed in source function (for example)
                dxdxi = 0.0; % Jacobian (matrix) at this quadrature point
                for ishl = 1:nshl
                  % xglobal = xglobal+Nq(ishl)*xl(ishl);
                  dxdxi = dxdxi+dNq(ishl)*xl(ishl);
                end
                detJ = dxdxi; % for 1D mesh
                WdetJ = wq(iq)*detJ;
                 % detJ is h_mesh/2
                hhalf_mesh = detJ;
                dxidx = 1.0/dxdxi; % Jacobian (matrix) inverse for 1D mesh
                for ishl = 1:nshl
                  dNdxq(ishl) = dNq(ishl)*dxidx; % global (or x) derivative of shape function at this quadrature point
                end

                % prob. parameters at this quadrature point
                % stab. parameters
                kappaq = kappa;
                if method==1
                    tauadv_invsq = (1/hhalf_mesh)^2*gPC_multiply(beta_cap,beta_cap,M,Ckij);
                    taudiff_invsq = (1/hhalf_mesh^2)^2*(kappaq)^2*const;
                    tau_invsq = tauadv_invsq + 9*taudiff_invsq;
                    tau_sq = gPC_divide(const,tau_invsq,M,Ckij);
                    tauq_cap = gPC_squareroot(tau_sq,M,Ckij,0,iel);  
                    tauq_cap = tauq_cap';
                else
                    tauq_cap = 0*const;
                end 

                sourceq = 1; % source term
                
                % Expectation calculation
                beta_tau_cap = gPC_multiply(beta_cap,tauq_cap,M,Ckij);
                beta_tau_beta_cap = gPC_multiply(beta_cap,beta_tau_cap,M,Ckij); 
                kapp_cap = kappaq*const;                
                exp_phiI = Ckij(1,I,1);
                exp_phiIJ = Ckij(1,I,J);
                exp_beta_tau_beta_phiIJ = 0;
                exp_beta_phiIJ = 0;
                exp_kapp_phiIJ = 0;
                for K=1:M
                    exp_kapp_phiIJ = exp_kapp_phiIJ + Ckij(K,I,J)/(2*K-1)*kapp_cap(K);
                    exp_beta_phiIJ = exp_beta_phiIJ + Ckij(K,I,J)/(2*K-1)*beta_cap(K);
                    exp_beta_tau_beta_phiIJ = exp_beta_tau_beta_phiIJ + Ckij(K,I,J)/(2*K-1)*beta_tau_beta_cap(K);                   
                end
                exp_beta_tau_phiI = 0;
                for K1=1:M
                    for K2=1:M
                        exp_beta_tau_phiI = exp_beta_tau_phiI + Ckij(I,K1,K2)/(2*I-1)*beta_cap(K1)*tauq_cap(K2);
                    end
                end
                
                %%%% prep data - end

                for ishl = 1:nshl
                  % VMS
                  be(ishl) = be(ishl)+((Nq(ishl)*exp_phiI + (dNdxq(ishl)*exp_beta_tau_phiI))*sourceq)*WdetJ;
                  for jshl = 1:nshl
                    % VMS 
                    Ae(ishl,jshl) = Ae(ishl,jshl)+(Nq(ishl)*exp_beta_phiIJ*dNdxq(jshl) + ...
                    dNdxq(ishl)*kappaq*dNdxq(jshl)*exp_phiIJ + ...
                    (dNdxq(ishl)*exp_beta_tau_beta_phiIJ*dNdxq(jshl)))*WdetJ;
                  end
                end
              %============
              end % end of loop over quadrature points
              
              % assemble Ae to A and be to b
              for ishl = 1:nshl
                for jshl = 1:nshl
                  A(ien(iel,ishl),ien(iel,jshl)) = A(ien(iel,ishl),ien(iel,jshl))+Ae(ishl,jshl);
                end
                b(ien(iel,ishl)) = b(ien(iel,ishl))+be(ishl);
              end

            %==============
            end % end of loop over elements
            K_cell{I,J} = A;
            F_cell{I,1} = b;
        end
    end
    % Boundary conditions
    for I=1:M
        btmp=F_cell{I,1};
        btmp(1)=g1*Ckij(1,I,1);
        btmp(npoints)=g2*Ckij(1,I,1);
        for J=1:M
            btmp(2)=btmp(2)-K_cell{I,J}(2,1)*g1*Ckij(1,J,1);
            btmp(npoints-1)=btmp(npoints-1)-K_cell{I,J}(npoints-1,npoints)*g2*Ckij(1,J,1);
        end
        F_cell{I,1}=btmp;
    end
    for I=1:M
        for J=1:M
            Atmp = K_cell{I,J};
            if I==J
                Atmp(1,1)=1;
                Atmp(npoints,npoints)=1;
                Atmp(1,2)=0;
                Atmp(2,1)=0;
                Atmp(npoints-1,npoints)=0;
                Atmp(npoints,npoints-1)=0;
            else
                Atmp(1,:)=0;
                Atmp(:,1)=0;
                Atmp(npoints,:)=0;
                Atmp(:,npoints)=0;
            end
            K_cell{I,J}=Atmp;
        end
    end
    
    % Final assmbly to matrix and vector
    K_matrix = cell2mat(K_cell);
    F_vector = cell2mat(F_cell);
    U_vector = K_matrix\F_vector;
    
    u = zeros(M,npoints);
    k=0;
    for j=1:M
        for i=1:npoints
            k=k+1;
            u(j,i) = U_vector(k);
        end
    end 
    cellpec = 2*hhalf_mesh/kappa; % beta_max = 2
    y=linspace(-1,1,y_points);
    legP = zeros(M,y_points);
    for J=1:M
        legP(J,:) = legendreP(J-1,y);
    end
    soln = zeros(npoints,y_points);
    for i=1:npoints
        for j=1:y_points
            for J=1:M
                soln(i,j) = soln(i,j) + u(J,i)*legP(J,j);
            end
        end            
    end
    figure(method);
    surface(xmesh,y,soln')
%     colormap gray
    axis([0 1 0 1 0 3])
    xlabel('x')
    ylabel('y')
    title(['P=',num2str(P),' \kappa=',num2str(kappa),' and N_{el}=',num2str(nel),' CellPec_{max}=',num2str(cellpec)])
    if method==1
        filename1=['Stoch_AdvDiff',num2str(q),'a'];
    else
        filename1=['Stoch_AdvDiff',num2str(q),'b'];
    end
%     print(filename1, '-dpng', '-r600');

    % expectation of u(x,y)
    expec = zeros(1,npoints);
    for i=1:npoints
        for J=1:M
            expec(i) = expec(i)+u(J,i)*Ckij(1,1,J);
        end
    end
    figure(3);
    if method==1
        plot(xmesh,expec,'k-','Linewidth',1.5)
        hold on
    else
        plot(xmesh,expec,'-.r','Linewidth',1)
    end
    
    % variance of u(x,y)
    var = 0*expec;
    for i=1:npoints
        for I=1:M
            for J=1:M
                var(i) = var(i) + u(I,i)*u(J,i)*(Ckij(1,I,J)-Ckij(1,I,1)*Ckij(1,J,1));
            end
        end
    end
    figure(4);
    if method==1
        plot(xmesh,var,'k-','Linewidth',1.5)
        hold on
    else
        plot(xmesh,var,'-.r','Linewidth',1)
    end    
end
figure(3)
title('Expectation')
legend('VMS','Galerkin','Location','northwest')
filename2 = ['Expectation',num2str(q)];
% print(filename2, '-dpng', '-r600');
hold off

figure(4)
title('Variance')
legend('VMS','Galerkin','Location','northwest')
filename3 = ['Variance',num2str(q)];
% print(filename3, '-dpng', '-r600');
hold off



    







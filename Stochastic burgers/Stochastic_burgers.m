clear;
clc;

% problem parameters
P = 4; % maximum gPC order
d = 1; % number of stochastic variable
M = factorial(P+d)/factorial(P)/factorial(d); %total number of terms in Pascal's simplex
mu = 0.05;
nr_limit=100;
[indx,Ckij] = lexico_table(P,d,1);
const = zeros(1,M);
const(1) = 1;% gPC coefficinets for const=1 [1 0 0 ... 0]

q=1;
% Stochastic BCs
g1 = 1*const; % deterministic 
g2 = g1; % stochastic   
g2(1) = 0;
g2(2) = 1;

source = 1;
nel = 100;
y_points = 50;
yy=linspace(-1,1,y_points);
nrlimit=100;

% num. of points
npoints = nel+1;

% domain (starting at x_A and of length L)
x_A = -1.0;
L = 2.0;
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

 
K_cell =  cell(M,M); % Jacobian
F_cell = cell(M,1); % Residual

for method=1:1
    % trial solution
    
    Utrial = cell(M,1);
    updt = cell(M,1);
    for I=1:M
        Utrial{I,1}=ones(npoints,1);
%         Utrial{I,1}=guess(:,I);
        updt{I,1}=ones(npoints,1);
        Utrial{I,1}(1)=g1(I);
        Utrial{I,1}(npoints)=g2(I);
    end
    nr=0; %newton-raphson iteration number
    updt_vector = cell2mat(updt);
    
    while (max(abs(updt_vector))>(1E-3))
        disp(max(abs(updt_vector)));
        nr = nr+1;
            if nr>nrlimit
                break;
            end
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
                  
                  % element level U_1^K and U_2^K
                  Unode = cell(M,1);
                  for K=1:M
                    Unode{K,1} = [Utrial{K,1}(ien(iel,1)),Utrial{K,1}(ien(iel,2))];
                  end

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
                    
                    alphaq_cap = 0*const;
                    betaq_cap = 0*const;
                    for K=1:M
                        for ishl=1:nshl
                            alphaq_cap(1,K) = alphaq_cap(1,K) + Unode{K,1}(ishl)*Nq(ishl);
                            betaq_cap(1,K) = betaq_cap(1,K) + Unode{K,1}(ishl)*dNdxq(ishl);
                        end
                    end

                    % prob. parameters at this quadrature point
                    % stab. parameters
                    kappaq = mu;
                    if method==1
                        tauadv_invsq = (1/hhalf_mesh)^2*gPC_multiply(alphaq_cap,alphaq_cap,M,Ckij);
                        taudiff_invsq = (1/hhalf_mesh^2)^2*(kappaq)^2*const;
                        tau_invsq = tauadv_invsq + 9*taudiff_invsq;
                        tau_sq = gPC_divide(const,tau_invsq,M,Ckij);
                        tauq_cap = gPC_squareroot(tau_sq,M,Ckij,0,iel);  
                        tauq_cap = tauq_cap';
                    else
                        tauq_cap = 0*const;
                    end 
                    % gPC approximations
                    thetaq_cap = gPC_multiply(alphaq_cap,tauq_cap,M,Ckij); % U*tau
                    gammaq_cap = gPC_multiply(thetaq_cap,betaq_cap,M,Ckij); % U*tau*dUdx
                    deltaq_cap = gPC_multiply(thetaq_cap,alphaq_cap,M,Ckij); % U*tau*U

                    sourceq = source; % source term

                    % Expectation calculation
                    exp_IJ = Ckij(1,I,J);
                    exp_I = Ckij(1,I,1);
                    alpha=0;beta=0;betaU=0;gamma=0;tau=0;
                    delta=0;alphabeta=0;betadelta=0;alphatau=0;
                    for K=1:M
                        alpha = alpha + alphaq_cap(K)*Ckij(K,I,J)/(2*K-1);
                        beta = beta + betaq_cap(K)*Ckij(K,I,J)/(2*K-1);
                        gamma = gamma + gammaq_cap(K)*Ckij(K,I,J)/(2*K-1);
                        delta = delta + deltaq_cap(K)*Ckij(K,I,J)/(2*K-1);
                        tau = tau + tauq_cap(K)*Ckij(K,I,J)/(2*K-1);
                        betaU = betaU + betaq_cap(K)*Ckij(K,I,1)/(2*K-1);
                        for H=1:M
                            alphabeta = alphabeta + alphaq_cap(K)*betaq_cap(H)*Ckij(K,I,H)/(2*K-1);
                            betadelta = betadelta + betaq_cap(K)*deltaq_cap(H)*Ckij(K,I,H)/(2*K-1);
                            alphatau = alphatau + alphaq_cap(K)*tauq_cap(H)*Ckij(K,I,H)/(2*K-1);
                        end
                    end
                        

                    %%%% prep data - end
                    for ishl = 1:nshl
                        be(ishl) = be(ishl) + (dNdxq(ishl)*mu*betaU + Nq(ishl)*alphabeta +...
                            dNdxq(ishl)*betadelta - Nq(ishl)*sourceq*exp_I - Nq(ishl)*sourceq*alphatau)*WdetJ;
                        for jshl = 1:nshl                 
                            Ae(ishl,jshl) = Ae(ishl,jshl) + (dNdxq(ishl)*mu*dNdxq(jshl)*exp_IJ...
                               + Nq(ishl)*Nq(jshl)*beta + Nq(ishl)*dNdxq(jshl)*alpha + ...
                               Nq(jshl)*dNdxq(ishl)*gamma + dNdxq(ishl)*Nq(jshl)*gamma + ...
                               dNdxq(ishl)*dNdxq(jshl)*delta - Nq(jshl)*dNdxq(ishl)*sourceq*tau)*WdetJ;
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
            F_cell{I}(1)=0;
            F_cell{I}(npoints)=0;
        end


        % Final assmbly to matrix and vector
        K_matrix = cell2mat(K_cell);
        F_vector = cell2mat(F_cell);
        updt_vector = K_matrix\F_vector;
        U_vector = cell2mat(Utrial);
        U_vector = U_vector-updt_vector;
        dim = npoints*ones(M,1);
        Utrial = mat2cell(U_vector,dim,1);        
    end
    
    
    u = zeros(M,npoints);
    k=0;
    for j=1:M
        for i=1:npoints
            k=k+1;
            u(j,i) = U_vector(k);
        end
    end 
    
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
    surface(xmesh,y,soln')%,'Edgecolor','none'
%     colormap gray
    view(30,30)
    xlabel('x')
    ylabel('y')
    title(['P=',num2str(P),' \mu=',num2str(mu),' and N_{el}=',num2str(nel)])
    if method==1
        filename1=['Stoch_Burgers',num2str(q),'a'];
    else
        filename1=['Stoch_Burgers',num2str(q),'b'];
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



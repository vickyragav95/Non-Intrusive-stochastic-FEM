%generates index of each univariate basis's power based on lexicographic
%ordering
%input:             p: max gPC order
%                   d: number of stochastic variables/dimensions
%output:            indx: matrix where each row is an entry of Pascal's simplex in lexicographic ordering
%                   Ckij: triple product tensor integral


function [indx,Ckij]=lexico_table(p,d,Ckij_flag)

M=factorial(p+d)/factorial(p)/factorial(d); %total number of terms in Pascal's simplex

indx=zeros(M,d); %matrix is M number of terms x d number of dimensions

%p=0 level
indx(1,1:d)=zeros(1,d); %indexing starts at 0, so 1st level is zero vector

%p=1 level is just p=0 level + unit vector for each d
M_prev=factorial(1+d)/factorial(d);
num_row=2;
for i=1:d
    e_i=unit_vector(d,i,0);
    indx(num_row,1:d)=e_i;
    num_row=num_row+1;
end
terms=ones(1,d); %used below
num_terms_prev=d; %used below, number of terms for p=1 level is number of dimensions
start_row=2; %starting row in indx for this level

%p>=2 levels
for i=2:p
    M_current=factorial(i+d)/factorial(i)/factorial(d); %total number of terms through current level
    num_terms=M_current-M_prev; %number of additional terms from this level
    for j=1:d %add unit vector appropriately to entries from previous level to get entries for this level
        tmp=terms(j);
        terms(j)=num_terms_prev;
        for k=1:terms(j)
            e_i=unit_vector(d,j,0);
            indx(num_row,1:d)=indx(start_row+k-1,1:d)+e_i;
            num_row=num_row+1;
        end
        num_terms_prev=num_terms_prev-tmp;
        start_row=start_row+tmp;
    end
    num_terms_prev=num_terms; %store for next level
    start_row=M_prev+1; %starting row in indx for this level
    M_prev=M_current; %store for next level
end

%Ckij calculation
if Ckij_flag==1
    Ckij=zeros(M,M,M);
    for i=1:M
        for j=1:M
            for k=1:M
                tmp_Ckij=1.0;
                for dim=1:d
                    i_dim=indx(i,dim);
                    j_dim=indx(j,dim);
                    k_dim=indx(k,dim);
                    g=(i_dim+j_dim+k_dim)/2.0;
                    if mod(2*g,2)==0 && abs(i_dim-j_dim)<=k_dim && k_dim<=i_dim+j_dim
                        tmp1=factorial(2*g-2*i_dim)*factorial(2*g-2*j_dim)*factorial(2*g-2*k_dim)/factorial(2*g+1);
                        tmp2=factorial(g)/factorial(g-i_dim)/factorial(g-j_dim)/factorial(g-k_dim);
                        tmp3=tmp2*tmp2;
                        three_jm_sym_sq=tmp1*tmp3;
                    else
                        three_jm_sym_sq=0.0;
                    end
                    tmp_three_jm=three_jm_sym_sq*(2*k_dim+1);%normalize triple product integral by (2k+1)
                    if abs(tmp_three_jm)<1.0e-14
                        tmp_Ckij=0.0;
                        break
                    end
                    tmp_Ckij=tmp_Ckij*tmp_three_jm;
                end
                Ckij(k,i,j) = tmp_Ckij;
            end
        end
    end
end


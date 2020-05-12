%performs gPC square root

function c=gPC_squareroot(a,M,Ckij,output,iel)

  if length(a)~=M
      error('a  must be of size M');
  end
  
  c=zeros(1,M);
  K=zeros(M,M);
  f=zeros(M,1);
  %set initial guess to sqrt of 1st coefficient
  b=[sign(a(1)).*sqrt(abs(a(1))),zeros(1,M-1)];
  f_L2=0.0; %L2 norm of initial RHS
  for k=1:M
      for l=1:M
          for j=1:M
              f(k)=Ckij(k,j,l)*b(j)*b(l)+f(k);
          end
      end
      f(k)=-f(k)+a(k);
      f_L2=f(k)*f(k)+f_L2;
  end
  f_L2=sqrt(f_L2);
  if output==1
      %fprintf('step: 0, L2 RHS: %.16e\n', f_L2)
  end
  f_L2_vals=f_L2;
  delb_L2_vals=[];
  %tolerances
  abs_del=f_L2; %used for solution update
  abs_tol=1.0e-12; %tolerance for nonlinear solve
  rel_del=1;
  rel_tol=1.0e-4;
  %step settings
  n=1; %current step number
  n_max=1000; %max step number
  stop_flag=0;
  stop_counter=0;
  %solve the nonlinear system
  while stop_flag<1 %newtons method
      if n<=n_max
          for k=1:M
              f(k)=0.0;
              for l=1:M
                  K(k,l)=0.0;
                  for j=1:M
                      K(k,l)=Ckij(k,j,l)*b(j)+Ckij(k,l,j)*b(j)+K(k,l);
                      f(k)=Ckij(k,j,l)*b(j)*b(l)+f(k);
                  end
              end
              f(k)=-f(k)+a(k); %new RHS
          end
          delb=K\f; %solve system to get update
          %lin_res=K*delb-f;
          %lin_res_L2=0.0; %L2 norm of linear solver residual
          delb_L2=0.0; %L2 norm of update
          for k=1:M
              %lin_res_L2=lin_res(k)*lin_res(k)+lin_res_L2;
              delb_L2=delb(k)*delb(k)+delb_L2; %calculate L2 norm of updated guess
              b(k)=b(k)+delb(k); %update guess
          end
          %lin_res_L2=sqrt(lin_res_L2);
          delb_L2=sqrt(delb_L2);
          if output==1
              %fprintf('step: %d, L2 linear residual: %.16e\n', n, lin_res_L2);
              %fprintf('step: %d, L2 update: %.16e\n', n, delb_L2);
          end
          delb_L2_vals=[delb_L2_vals; delb_L2];
          f_L2_old=f_L2;
          f_L2=0.0; %L2 norm of RHS
          for k=1:M
              f(k)=0.0;
              for l=1:M
                  for j=1:M
                      f(k)=Ckij(k,j,l)*b(j)*b(l)+f(k);
                  end
              end
              f(k)=-f(k)+a(k);
              f_L2=f(k)*f(k)+f_L2; %calculate L2 norm of RHS
          end
          f_L2=sqrt(f_L2);
          if output==1
              %fprintf('step: %d, L2 updated RHS: %.16e\n', n, f_L2)
          end
          f_L2_vals=[f_L2_vals; f_L2];
          %res=a-gPC_multiply(b,b,M,Ckij);
          %rel_del=max(abs(res))/max(abs(a));
          %abs_del=max(abs(res));
          if f_L2/f_L2_old > 1.0/rel_tol %if spike, revert to previous and stop
              stop_counter=stop_counter+1;
              if stop_counter>1
                  b(k)=b(k)-delb(k); %undo update
                  f_L2_vals=f_L2_vals(1:end-1);
                  delb_L2_vals=delb_L2_vals(1:end-1);
                  break
              end
          end
          rel_del=f_L2/f_L2_old;
          %disp(rel_del);
          abs_del=f_L2;
          if rel_del<rel_tol || abs_del<abs_tol
              stop_flag=1;
          end
          n=n+1; %increment counter
      else
          break
      end
      c=b';
  end
  if output==1
    figure
    semilogy(0:n-1,f_L2_vals);
    set(gca,'FontSize',20);
    title(sprintf('f_L2 for element %d', iel));
    
    figure
    semilogy(1:n-1,delb_L2_vals);
    set(gca,'FontSize',20);
    title(sprintf('delb_L2 for element %d', iel));
  end  
  %disp(del);
%performs gPC truncated multiplication

function c=gPC_multiply(a,b,M,Ckij)

  if length(a)~=M || length(b)~=M
      error('a and b must be of size M');
  end
  
  c=zeros(1,M);
  for k=1:M
      for i=1:M
          for j=1:M
              c(k)=c(k)+Ckij(k,i,j)*a(i)*b(j);
          end
      end
  end
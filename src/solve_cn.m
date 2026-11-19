function [cn]=solve_cn(N)
    cn=zeros(1,N);
    for n=1:N
        k=((-1)^(n+1))/(2*n-1);
        cn(n)=(k*fuzhu1(n,N))/(fuzhu2(n)*fuzhu3(n,N));
    end
end
function [a1]=fuzhu1(n,N)
    a1=1;
    for l=1:N
        if l~=n
            a1=a1*(2*l-1)^2;
        end
    end
end
function [a2]=fuzhu2(n)
    a2=1;
    for l=1:n-1
        a2=a2*((2*n-1)^2-(2*l-1)^2);
    end
end
function [a3]=fuzhu3(n,N)
    a3=1;
    for l=n+1:N
        a3=a3*((2*l-1)^2-(2*n-1)^2);
    end
end

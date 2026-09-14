%Smolyak with symmetric and disjoint sets,
%Parisa Kamali
%Oct 14,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function[psi]=Psi(d,mu,y)

%initial condition
%-------------------------------------------------------------------------------
d=2;                            %dimension
mu=2;                           %approximation level
y=[1;0];
n=size(y,2);
%building the whole matrix V_T given d and \mu  
%-------------------------------------------------------------------------------
v_0=ones(d,1);                  %initial vector for i_1,...1_d
V_T=v_0;                        %V_T is the vector that will contain all combinations of i_1,..., i_d
if mu==0                        %in \mu==0 we are doen already
    break;
else                                                   
    k=1;
    v=[];
                        
    for l=1:d                          
       v_0(l,1)=v_0(l,1)+k;
       v=[v v_0];
       v_0(l,1)=v_0(l,1)-k;
    end
    
    v_T=v;
    V_T=[V_T v]; 
    k=k+1;

    while k<mu+1 
       v=v_T; 
       v_T=[];
       for l=1:size(v,2)
           for j=1:d
               v(j,l)=v(j,l)+1;
               if isempty(v_T)==0
                    t=ismember(v(:,l)',v_T','rows');
                    if t==0
                       v_T=[v_T v(:,l)];
                    end
               else
                   v_T=[v_T v(:,l)];
               end
               v(j,l)=v(j,l)-1;
           end
       end
       k=k+1;
       V_T=[V_T v_T];
    end
end 






%constructing vector A containing \psi
%-----------------------------------------------------------------------------------
% constructing the matrix that associates the number of \psi in each
% dimension
D=zeros(size(V_T,1),size(V_T,2));
P=zeros(1,size(V_T,2));

for j=1:size(V_T,2)
    for k=1:size(V_T,1)
        if V_T(k,j)==1
            D(k,j)=1;
        elseif V_T(k,j)==2
            D(k,j)=3;
        else
            D(k,j)=3;
            for l=3:V_T(k,j)
                D(k,j)=D(k,j)+2^(l-2);
            end
            
        end
    end
     P(1,j)=prod(D(:,j));
end
PSI=[];            %size of the final matrix

% number of \psi in each dimension
%--------------------------------------------------------------------------
m=max(V_T');                %will be used in construction of array of \psi
M=max(D');

 psi_T=zeros(d,max(M));
 
for j=1:max(M)
    if j==1
        psi_T(:,1)=1;
    elseif j==2
        psi_T(:,2)=y(:);
    else 
        psi_T(:,j)=2.*y.*psi_T(:,j-1)-psi_T(:,j-2);
    end
end
        
for j=1:size(V_T,2)      
    psi=Inf(d,max(2,2^(max(m)-2)));
    for k=1:d
        if V_T(k,j)==1
            psi(k,1)=1;
        elseif V_T(k,j)==2
            psi(k,1:2)=[y(k) 2*y(k)*y(k)-1];
        else
            b=psi_T(k,2^(V_T(k,j)-2)+2:2^(V_T(k,j)-2)+2+2^(V_T(k,j)-2)-1);
            psi(k,1:length(b))=b;
            
        end
    end
    a=1;
    for k=1:d
        a=kron(a,psi(k,:));
    end
    for k=length(a):-1:1
        if isnan(a(k))==1 || isinf(a(k))==1
            a(k)=[];
        end
    end
    PSI=[PSI;a'];
end

        









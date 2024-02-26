function [IGB,X,Y]=Ince_Gaussian(L,N,parity,p,m,e,w0,k,z)
%{
Calculate an Ince-Gaussian Beam at a given z plane

INPUTS:
    L = transverse physical size of the X-Y space [-L,L,-L,L]
    N = number of sampling points (must be ODD)
    parity = parity of the beam, 0 = EVEN  C;  1 = ODD  S
    p, m = order and degree of the Ince Gaussian beam
           p=0,1,2,3,... for parity=0 and p=1,2,3,.. for parity=1
           0<=m<=p   
           (p,m) must have the same parity, i.e., (-1)^(p-m)=1
    e = ellipticity parameter
    w0 = beam width(waist) at z=0
    k = 2*pi/lambda, wavenumber, lambda=wavelength 
    z = propagation distance

OUTPUTS:
    IGB = Ince Gaussian Beam (Intensity is analytically normalized to 1, i.e., Integral(|IGB|^2)=1)
    X and Y = space matrices

EXAMPLE:   
        [IGB,X,Y]=Ince_Gaussian(15e-3,501,0,6,2,2,3e-3,(2*pi/632.8e-9),0); figure;
        surf(X,Y,abs(IGB)); shading interp; lighting phong; view(2); axis tight; axis equal; xlabel('x','Fontsize',12); ylabel('y','Fontsize',12); title('Ince Gausssian Beam','Fontsize',12)
        dx=abs(X(2,2)-X(2,1)); Normalization=sum(sum(IGB.*conj(IGB))).*dx^2
    
For more information about Ince Gaussian Beams:

    "Ince-Gaussian beams," Miguel A. Bandres and J. C. Gutierrez-Vega
    Optics Letters, 29(2), 144-146 (2004) ( http://goo.gl/U18mol )

    "Ince-Gaussian modes of the paraxial wave equation and stable resonators," Miguel A. Bandres and Julio C. Gutierez-Vega
    Journal of the Optical Society of America A, 21(5), 873-880 (2004) ( http://goo.gl/rqq7nQ )

    "Observation of Ince-Gaussian modes in stable resonators," Ulrich T. Schwarz, Miguel A. Bandres and Julio C. Gutierrez-Vega
    Optics Letters, 29(16), 1870-1872 (2004) ( http://goo.gl/7lkSb4 )

AUTHOR: Miguel A. Bandres
URL:    www.mabandres.com
%}

%% CHECK INPUT
if mod(N,2)==0;      error('ERROR: N must be ODD'); end;
if parity==0
    if (m<0)||(m>p); error('ERROR: Wrong range for "m", 0<=m<=p'); end;
else
    if (m<1)||(m>p); error('ERROR: Wrong range for "m", 1<=m<=p'); end;  
end;
if (-1)^(m-p)~=1;    error('ERROR: (p,m) must have the same parity, i.e., (-1)^(m-p)=1');  end;

%% PARAMETERS
f0 = sqrt(e/2)*w0; % Focal distance of elliptic coordinates at z=0

%% INCE GAUSSIAN BEAM
if z==0 % At the z=0
   [xhi,etha,X,Y]=mesh_elliptic(f0,L,N); % Calculate elliptic coordinates (xhi,etha), (X,Y)= Cartesian Coordinates
   R = sqrt(X.^2 + Y.^2);                % Radial coordinate
   if parity==0
      IGB=CInceIGB(p,m,e,etha).*CInceIGB(p,m,e,1i*xhi).*exp(-(R/w0).^2);  % Even IGB
   else
      IGB=SInceIGB(p,m,e,etha).*SInceIGB(p,m,e,1i*xhi).*exp(-(R/w0).^2);  % Odd IGB
   end
else % At z~=0
   zr=1/2*k*w0^2;            % Rayleigh range 
   wz=w0*sqrt(1+(z./zr).^2);  % Beam width(waist) at z=0   
   Rz=z.*(1+(zr./z).^2);      % Radius of curvature of the phase front
   f = f0*wz/w0;             % Focal distance of elliptic coordinates at z   
   [xhi,etha,X,Y]=mesh_elliptic(f,L,N); % Calculate elliptic coordinates (xhi,etha), (X,Y)= Cartesian Coordinates
   R = sqrt(X.^2 + Y.^2);    % Radial coordinate
   if parity==0
      IGB=(w0/wz)*(CInceIGB(p,m,e,etha).*CInceIGB(p,m,e,1i*xhi)).*exp(-(R/wz).^2).* ...
         exp(1i*(k*z + k*R.^2/(2*Rz)-(p+1)*atan(z/zr)));    % Even IGB
   else
      IGB=(w0/wz)*(SInceIGB(p,m,e,etha).*SInceIGB(p,m,e,1i*xhi)).*exp(-(R/wz).^2).* ...
         exp(1i*(k*z + k*R.^2/(2*Rz)-(p+1)*atan(z/zr)));    % Odd IGB
   end   
end

%% COMPUTE THE NORMALIZATION CONSTANTS
if parity == 0;  
   if mod(p,2)==0;   
      [C0,~,coef]=CInceIGB(p,m,e,0);
      [Cp,~,~]=CInceIGB(p,m,e,3.14159/2);
      Norm = (-1)^(m/2)*sqrt(2)*gamma(p/2+1)*coef(1) *sqrt(2/3.14159)/w0/C0/Cp;  % Calculate normalization constant
   else
      [C0,~,coef]=CInceIGB(p,m,e,0); 
      [~,~,~,DCp]=CInceIGB(p,m,e,3.14159/2);
      Norm = (-1)^((m+1)/2) * gamma((p+1)/2+1) * sqrt(4*e/3.14159) * coef(1) / w0 / C0 / DCp; % Calculate normalization constant
   end
else             
   if mod(p,2)==0;   
      [~,~,coef,dS0]=SInceIGB(p,m,e,0);
      [~,~,~,dSp]=SInceIGB(p,m,e,3.14159/2); 
      Norm = (-1)^(m/2)*sqrt(2)*e*gamma((p+2)/2+1)*coef(1) *sqrt(2/3.14159)/w0/ dS0 / dSp;  % Calculate normalization constant      
   else
      [Sp,~,coef,~]=SInceIGB(p,m,e,3.14159/2); 
      [~,~,~,dS0]=SInceIGB(p,m,e,0);       
      Norm = (-1)^((m-1)/2) * gamma((p+1)/2+1) * sqrt(4*e/3.14159) * coef(1) / w0 / Sp / dS0;  % Calculate normalization constant      
   end
end

IGB = IGB*Norm; % Normalized the IGB to Integral(|IGB|^2)=1

return


%% SUBFUNCTIONS

%% Elliptic Coordinates
function [xi,eta,X,Y]=mesh_elliptic(f,L,N)
%{
Creates an elliptic coordinate mesh of size 2Lx2L with N points

INPUT:
    f = focal distance of elliptical coordinates 
    L = spatial size of the mesh, i.e., [-L,L]x[-L,L]
    N = ODD number of discrete points in the mesh. Must be ODD.

OUPUT:
   xi  = elliptic coordinate mesh 
   eta = eta elliptic coordinate mesh
   X   = x-Cartesian coordinate mesh
   Y   = y-Cartesian coordinate mesh
%}

%% Cartesian Coordinates
[X,Y]=meshgrid(linspace(-L,L,N));

%% Elliptic Coordinates
xi=zeros(N); eta=zeros(N);  % Initialization

% Calculate First Quadrant             
en = acosh((X(1:(N+1)/2,(N+1)/2:N)+1i*Y(1:(N+1)/2,(N+1)/2:N))/f);
ee = real(en);
nn = imag(en);
nn = nn + (nn<0)*2*3.14159; 
xi(1:(N+1)/2,(N+1)/2:N)=ee;
eta(1:(N+1)/2,(N+1)/2:N)=nn;
      
% Calculate other quadrants by symmetry       
xi(1:(N+1)/2 , 1:(N-1)/2)=fliplr(xi(1:(N+1)/2,(N+3)/2:N));
xi((N+3)/2:N,1:N)=flipud(xi(1:(N-1)/2,1:N));
eta(1:(N+1)/2,1:(N-1)/2)=3.14159-fliplr(eta(1:(N+1)/2,(N+3)/2:N));
eta((N+3)/2:N,1:N)=3.14159+rot90((eta(1:(N-1)/2,1:N)),2);

return

%% EVEN Ince Polynomial
function [IP,eta,coef,dIP]=CInceIGB(p,m,q,z) 
%{
This function calculates the EVEN Ince polynomials which are solutions of the Ince equation.
The Ince equation is a periodic linear second-order differential equation that has two families of independent solutions, 
namely, the even C^{m}_p(z,q) and odd S^{m}_p(z,q) Ince polynomials of order p and degree m. 
Physical considerations are such that Ince polynomials are periodic with period 2p. The values of \etha that satisfy this condition
are the eigenvalues of the Ince equation which is given by

d^2F/dz^2 + q*sin(2z)dF/dz+(\eta-p*q*cos(2z))*F=0

where
0<=z<2pi
p=0,1,2,3,...  
q=complex parameter
\etha = eigenvalue of the Ince equation.

INPUTS:
    p=0,1,2,3... Order of Ince polynomial
    0<=m<=p      m is the degree of the Ince polynomial
    (p,m) must have the same parity, i.e., (-1)^(p-m)=1
    q = complex parameter
    0<=z<2pi     independent variable (Vector or Matrix)

OUTPUTS:
    IP=C^{m}_p(z,q) EVEN Ince Polynomial
    eta  = eigenvalue of the Ince Polynomial 
    coef = coefficients of the Ince Polynomial
    dIP  = derivative of the Ince Polynomial

For more information about Ince Polynomial:

    Miguel A. Bandres and Julio C. Gutiérrez-Vega, "Ince-Gaussian modes of the paraxial wave equation and stable resonators"
    Journal of the Optical Society of America A, 21(5), 873-880 (2004) ( http://goo.gl/rqq7nQ )

AUTHOR: Miguel A. Bandres
URL:    www.mabandres.com
%}

%% Check Input
if (m<0)||(m>p); error('ERROR: Wrong range for "m", 0<=m<=p'); end;
if (-1)^(m-p)~=1;   error('ERROR: (p,m) must have the same parity, i.e., (-1)^(m-p)=1');  end;
[largo,ancho]=size(z); % change input to vector format
z=transpose(z(:));
normalization=1;

%% Calculate the Coefficients 
if mod(p,2)==0
    %%%% p Even %%%%
    j=p/2;  N=j+1;  n=m/2+1; 
    
    % Matrix
    M=diag(q*(j+(1:N-1)),1) + diag([2*q*j,q*(j-(1:N-2))],-1) + diag([0,4*((0:N-2)+1).^2]);
    if p==0; M=0; end;
        
    % Eigenvalues and Eigenvectors 
    [A,ets]=eig(M);
    ets=diag(ets); 
    [ets,index]=sort(ets);
    A=A(:,index);
    
    % Normalization
    if normalization==0;  
       N2=2*A(1,n).^2+sum(A(2:N,n).^2);
       NS=sign(sum(A(:,n)));
       A=A/sqrt(N2)*NS;
    else 
       mv=(2:2:p).';
       N2=sqrt(A(1,n)^2*2*gamma(p/2+1)^2+sum((sqrt(gamma((p+mv)/2+1).*gamma((p-mv)/2+1) ).*A(2:p/2+1,n)).^2 ));
       NS=sign(sum(A(:,n)));
       A=A/N2*NS;
    end
    
    % Ince Polynomial
    r=0:N-1;
    [R,X]=meshgrid(r,z);
    IP=cos(2*X.*R)*A(:,n);
    dIP=-2*R.*sin(2*X.*R)*A(:,n);
    eta=ets(n);
        
else
    %%%% p ODD %%%
    j=(p-1)/2;  N=j+1;  n=(m+1)/2; 
    
    % Matrix
    M=diag(q/2*(p+(2*(0:N-2)+3)),1)+diag(q/2*(p-(2*(1:N-1)-1)),-1) + diag([q/2+p*q/2+1,(2*(1:N-1)+1).^2]);

    % Eigenvalues and Eigenvectors 
    [A,ets]=eig(M);
    ets=diag(ets);
    [ets,index]=sort(ets);
    A=A(:,index);

    % Normalization
    if normalization==0;  
       N2=sum(A(:,n).^2);
       NS=sign(sum(A(:,n)));
       A=A/sqrt(N2)*NS;
    else
        mv=(1:2:p).';
        N2=sqrt(sum( ( sqrt(gamma((p+mv)/2+1).*gamma((p-mv)/2+1) ).*A(:,n)).^2 ));
        NS=sign(sum(A(:,n)));
        A=A/N2*NS;
    end
    
    % Ince Polynomial
    r=2*(0:N-1)+1;  
    [R,X]=meshgrid(r,z);
    IP=cos(X.*R)*A(:,n);
    dIP=-R.*sin(X.*R)*A(:,n);
    eta=ets(n);
end

coef=A(:,n);
IP=reshape(IP,[largo,ancho]); % reshape output to original format
dIP=reshape(dIP,[largo,ancho]);

return

%% ODD Ince Polynomial
function [IP,eta,coef,dIP]=SInceIGB(p,m,q,z)
%{
This function calculates the ODD Ince polynomials which are solutions of the Ince equation.
The Ince equation is a periodic linear second-order differential equation that has two families of independent solutions, 
namely, the even C^{m}_p(z,q) and odd S^{m}_p(z,q) Ince polynomials of order p and degree m. 
Physical considerations are such that Ince polynomials are periodic with period 2p. The values of \etha that satisfy this condition
are the eigenvalues of the Ince equation which is given by

d^2F/dz^2 + q*sin(2z)dF/dz+(\eta-p*q*cos(2z))*F=0

where
0<=z<2pi
p=0,1,2,3,...  
q=complex parameter
\etha = eigenvalue of the Ince equation.

INPUTS:
    p=0,1,2,3... Order of Ince polynomial
    0<=m<=p      m is the degree of the Ince polynomial
    (p,m) must have the same parity, i.e., (-1)^(p-m)=1
    q = complex parameter
    0<=z<2pi     independent variable (Vector or Matrix)

OUTPUTS:
    IP=S^{m}_p(z,q) ODD Ince Polynomial
    eta  = eigenvalue of the Ince Polynomial 
    coef = coefficients of the Ince Polynomial
    dIP  = derivative of the Ince Polynomial

For more information about Ince Polynomial:

    Miguel A. Bandres and Julio C. Gutiérrez-Vega, "Ince-Gaussian modes of the paraxial wave equation and stable resonators"
    Journal of the Optical Society of America A, 21(5), 873-880 (2004) ( http://goo.gl/rqq7nQ )

AUTHOR: Miguel A. Bandres
URL:    www.mabandres.com
%}

%% Check Input
if (m<1)||(m>p);  error('ERROR: Wrong range for "m", 1<=m<=p'); end;    
if (-1)^(m-p)~=1;   error('ERROR: (p,m) must have the same parity, i.e., (-1)^(m-p)=1');  end;
[largo,ancho]=size(z);  % change input to vector format
z=transpose(z(:));
normalization=1;

%% Calculate the Coefficients 
if mod(p,2)==0
    %%%% p Even %%%%
    j=p/2;  N=j+1; n=m/2; 
    
    % Matrix 
    M=diag(q*(j+(2:N-1)),1)+diag(q*(j-(1:N-2)),-1) + diag(4*((0:N-2)+1).^2);

    % Eigenvalues and Eigenvectors 
    [A,ets]=eig(M);
    ets=diag(ets);
    [ets,index]=sort(ets);
    A=A(:,index);

    % Normalization
    r=1:N-1;
    if normalization==0;  
       N2=sum(A(:,n).^2);
       NS=sign(sum(r.*transpose(A(:,n))));
       A=A/sqrt(N2)*NS;
    else 
        mv=(2:2:p).';
        N2=sqrt(sum((sqrt(gamma((p+mv)/2+1).*gamma((p-mv)/2+1) ).*A(:,n)).^2 ));
        NS=sign(sum(r.*A(:,n)'));        
        A=A/N2*NS;
    end
    
    % Ince Polynomial
    [R,X]=meshgrid(r,z);
    IP=sin(2*X.*R)*A(:,n);
    dIP=2*R.*cos(2*X.*R)*A(:,n);
    eta=ets(n);
    
else
    %%%% p ODD %%%
    j=(p-1)/2; N=j+1; n=(m+1)/2; 
    
    % Matrix
    M=diag(q/2*(p+(2*(0:N-2)+3)),1)+diag(q/2*(p-(2*(1:N-1)-1)),-1) + diag([-q/2-p*q/2+1,(2*(1:N-1)+1).^2]);

    % Eigenvalues and Eigenvectors  
    [A,ets]=eig(M);
    ets=diag(ets);
    [ets,index]=sort(ets);
    A=A(:,index);
    
    % Normalization
    r=2*(0:N-1)+1;  
    if normalization==0;  
       N2=sum(A(:,n).^2);
       NS=sign(sum(r.*transpose(A(:,n))));
       A=A/sqrt(N2)*NS;
    else
       mv=(1:2:p).';
       N2=sqrt(sum( ( sqrt(gamma((p+mv)/2+1).*gamma((p-mv)/2+1) ).*A(:,n)).^2 ));
       NS=sign(sum(r.*A(:,n)'));
       A=A/N2*NS;
    end
    
    % Ince Polynomial
    [R,X]=meshgrid(r,z);
    IP=sin(X.*R)*A(:,n);
    dIP=R.*cos(X.*R)*A(:,n);
    eta=ets(n);
        
end

coef=A(:,n);
IP=reshape(IP,[largo,ancho]);  % reshape output to original format
dIP=reshape(dIP,[largo,ancho]);

return





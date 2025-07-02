
!Programa para resolver la función de onda en un pozo cuadrado

program cuantica
implicit none

100 format(i10,3x,e25.17)

!Establecemos las constantes del programa
integer*8,parameter::N=1000,T_max=5000
real*8,parameter::n_ciclos=20.d0,s1=1.d0,h=0.5d0,lambda=1.d0,k_0=2.d0*dacos(-1.d0)*n_ciclos/dfloat(N),s=s1/(4.d0*k_0**2.d0)

!2000: 1700
!1000: 460
!500: 110


!n_c:10->100
!n_c:20->190
!n_c:50->500


!Parametros
real*8 norm,suma,V,pos,desv,dran_u
integer*8 j,k,l

!Constantes
complex*16,parameter::i=(0.d0,1.d0)

!Vectores en variable compleja
complex*16 phi(0:N),A_menos,A_mas,A_0(0:N),alpha(0:N),b(0:N),beta(0:N),chi(0:N)

!Inicializamos los vectores, aun que los de la posicion N no se van a emplear
do j=0,N
    phi(j)=complex(0.d0,0.d0)
    A_0(j)=complex(0.d0,0.d0)
    alpha(j)=complex(0.d0,0.d0)
    b(j)=complex(0.d0,0.d0)
    beta(j)=complex(0.d0,0.d0)
    chi(j)=complex(0.d0,0.d0)
end do

call dran_ini(1231231)

norm=0.d0


!Generamos la función de onda y calculamos su norma
do j=1,N-1
    phi(j)=complex(dcos(k_0*dfloat(j)),dsin(k_0*dfloat(j)))*dexp(-8.d0*(dfloat(4*j-N)/(dfloat(N)))**2.d0)
    norm=norm+phi(j)*dconjg(phi(j))*h
end do

do j=0,N
    !phi(j)=phi(j)/dsqrt(norm)
end do

!=Calculamos los valores de A
A_menos=complex(1.d0,0.d0)
A_mas=complex(1.d0,0.d0)


do j=0,N
    
    !Calculamos previamente el valor del potencial en j
    if((dfloat(j).le.3.d0*dfloat(N)/5.d0).and.(dfloat(j).ge.2.d0*dfloat(N)/5.d0)) then
        V=lambda*k_0**2.d0
    !En caso contrario devolvemos el valor contrario
    else
            V=0.d0
    end if
    
    A_0(j)=complex(-2.d0-V,2.d0/s)
    
end do


alpha(N-1)=complex(0.d0,0.d0)

!Calculamos los alpha
do j=N-2,0,-1
    alpha(j)=-A_menos/(A_0(j+1)+A_mas*alpha(j+1))
end do

!Empezamos la simulación
open(10,file='funcion_onda.dat')
open(11,file='norma.dat')



do j=0,T_max
    
   !Guardamos la norma y la funcion de onda
    if(mod(j,10).eq.0) then
       do k=0,N
          write(10,*)k,' ',Real(phi(k)*dconjg(phi(k)))
       end do
       write(10,*)
       write(10,*)

    end if

    suma=0.d0
    do k=0,N
        suma=suma+phi(k)*dconjg(phi(k))*h
    end do 
    
    
    write(11,100)j,suma
    write(11,*)
    write(11,*)
    
    !Calculamos los b
    do k=0,N
        b(k)=4.d0*i*phi(k)/s
    end do
    
    !Condicion de contorno de beta
    beta(N-1)=complex(0.d0,0.d0)
    
    !Calculamos los beta
    do k=N-2,0,-1
        beta(k)=(b(k+1)-A_mas*beta(k+1))/(A_0(k+1)+A_mas*alpha(k+1))
    end do

    !Calculamos los valores de chi
    chi(N)=complex(0.d0,0.d0)
    chi(0)=complex(0.d0,0.d0)
    
    do k=1,N-1
        chi(k)=alpha(k-1)*chi(k-1)+beta(k-1)
    end do
    
    !Calculamos la nueva funcion de onda
    do k=0,N
        phi(k)=chi(k)-phi(k)
    end do

end do

close(10)
close(11)



stop
end

include 'dranxor.f'

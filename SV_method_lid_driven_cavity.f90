program   SV_method
!This program solves the Lid Driven cavity problem through an implementation of the Streamfunction Vorticity method
implicit none
!defining  a square fluid flow domain
real,parameter::lx=1.0
real,parameter::ly=1.0
!defining grid
integer,parameter::nx=100
integer,parameter::ny=100
!creating arrays
real,dimension(0:nx+1,0:ny+1)::omega,psi,u,v ,domega_dt ,x,y,vort_rhs_laplacian,vort_rhs_advective ,omeganp1,psi_old,V_net,omega_old!streamfunction -psi,voricity -omega
real::dx,dy,uwall,nu,dt,Tp,h,psi_residual,vorticity_residual,tol
integer::i,j,nsteps,time

dx=lx/nx
dy=ly/ny
uwall=8 !top wall speed in m/s
nu=0.01 !kinematic viscosity
dt=0.001
Tp=6.0
nsteps=int(Tp/dt)
h=dx

!intitial guess for stream function and vorticity are set to zero
psi(0:nx+1,0:ny+1)=0
omega(0:nx+1,0:ny+1)=0
domega_dt(0:nx+1,0:ny+1)=0

!setting boundary conditions
omega(0:nx+1,ny+1)=(-2*uwall)/(dy)
vort_rhs_laplacian(0:nx+1,0:ny+1)=0
vort_rhs_advective(0:nx+1,0:ny+1)=0


do i=0,nx+1
 do j=0,ny+1
x(i,j)=i*dx
y(i,j)=j*dy
end do
end do
open(1,file="x_y_u_v_Vtotal.dat",status='replace')
open(2,file="velocity.plt",status='replace')
open(3,file="x_y_Psi_omega.dat",status='replace')
open(4,file="Psi.plt",status='replace')



do time=1,nsteps
    psi_old=psi
    omega_old=omega
    do i=0,nx+1
        do j=0,ny+1
            if(((i.NE.0).AND.(i.NE.nx+1)).AND.((j.NE.0).AND.(j.NE.ny+1)))then   !checking for interior points

                vort_rhs_laplacian(i,j)=nu*(((omega(i+1,j)-(2*omega(i,j))+omega(i-1,j))/(dx**2)) + (((omega(i,j+1)-(2*omega(i,j))+omega(i,j-1))/(dy**2))))
                
                vort_rhs_advective=(((psi(i,j+1)-psi(i,j-1))/(2*dy))*((omega(i+1,j)-omega(i-1,j))/(2*dx)))-(((psi(i+1,j)-psi(i-1,j))/(2*dx))*((omega(i,j+1)-omega(i,j-1))/(2*dy)))

                omega(i,j)=omega(i,j)+(vort_rhs_laplacian(i,j)-vort_rhs_advective(i,j))*dt
            
                omeganp1(i,j)=omega(i,j)
            else if((j.EQ.0).AND.( i.NE.0) .AND.(i.NE.nx+1))then     !bottom B.C
                 omega(i,j)=2*(psi(i,j)-psi(i,j+1))/dy**2;

            else if( (j.EQ.ny+1).AND. (i.NE.0) .AND.(i.NE.nx+1))then    !top B.Cs
                omega(i,j)=2*(psi(i,j)-psi(i,j-1))/(dy**2)-(2*uwall/dy)

            else if((i.EQ.0).AND. (j.NE.0) .AND. (j.NE.ny+1))then   !left wall B.C
                omega(i,j)=2*(psi(i,j)-psi(i+1,j))/(dx**2)
            
            else if( (i.EQ.nx+1).AND.( j.NE.0).AND.(j.NE.ny+1))then  !right wall BC
                omega(i,j)=2*(psi(i,j)-psi(i-1,j))/(dx**2)
            
            else if((i.EQ.0).AND.(j.EQ.0))then   !left corner
                omega(1,j)=0

            else if((i.EQ.nx+1).AND.(j.EQ.0))then !right corner
                   omega(1,j)=0

            else if((i.EQ.0).AND.(j.EQ.ny+1))then    !top left corner
                    omega(i,j)=(omega(i+1,j))

            else if((i.EQ.nx+1).AND.(j.EQ.ny+1))then  !top right corner
                     omega(i,j)=(omega(i-1,j))
            end if
        
        end do 
        
    end do
    
!computing stream function
    do i=0,nx+1
       do j=0,ny+1
        if ((j.NE.0 ).AND.(j.NE.ny+1).AND.(i.NE.0) .AND.(i.NE.nx+1))then
       psi(i,j)=0.25*((h**2)*omega(i,j)+psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1));
       else
       end if
       u(i,j)=1.0
       v(i,j)=0.0
       end do
    end do


    do i=0,nx+1
        do j=0,ny+1
        if ((j.NE.0).AND. (j.NE.ny+1) .AND. (i.NE.0).and.(i.NE.nx+1))then
            u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*dy)
             v(i,j)=-(psi(i+1,j)-psi(i-1,j))/(2*dx)
            else if(j.EQ.ny+1)then
            u(i,j)=1;
            v(i,j)=0;
            else
                u(i,j)=0;
                v(i,j)=0;
        end if
        V_net(i,j)=SQRT(u(i,j)**2+v(i,j)**2)
        end do

    end do
 
psi_residual=MAXVAL(ABS(psi-psi_old))
vorticity_residual=MAXVAL(ABS(omeganp1-omega_old))


 tol=0.0001  !convergence criteria


write(*,*)
write(*,*)"The time step is:",time
write(*,*)
write(*,*)"Streamfunction residual:-"
write(*,*)psi_residual
write(*,*)
write(*,*)"vorticity residual:-"
write(*,*)vorticity_residual
write(*,*)
write(*,*)

if(psi_residual<tol)then
write(*,*)"convergence of",tol,"after",time,"steps,exiting solution"
exit
elseif(time.EQ.nsteps)then
write(*,*)"Max iterations reached,exiting solver......"
end if


end do
!writing file after  end of solution
do i=0,nx+1
do j=0,ny+1
write(1,*)x(i,j),y(i,j),u(i,j),v(i,j),V_net(i,j)
write(3,*)x(i,j),y(i,j),psi(i,j),omega(i,j)
end do
end do
write(2,*)'set xlabel "x"'
write(2,*)'set ylabel "y"'
write(2,*)'set zlabel "V"'
write(2,*)'set title "velocity"'
write(2,*)'set autoscale xy'
write(2,*)'plot "x_y_u_v_Vtotal.dat" using 1:2:5 with image'
CALL SYSTEM('gnuplot -p velocity.plt')

write(4,*)'set xlabel "x"'
write(4,*)'set ylabel "y"'
write(4,*)'set zlabel "Stream function"'
write(4,*)'set title "Psi(Streamfunction)"'
write(4,*)'set autoscale xy'
write(4,*)'plot "x_y_Psi_omega.dat" using 1:2:3 with image'
CALL SYSTEM('gnuplot -p Psi.plt')

close(1)
close(2)
close(3)

end program
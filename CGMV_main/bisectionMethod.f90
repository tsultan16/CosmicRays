!Bisection-Method Root Finder Algorithm

program bisectionRoot
implicit none

integer::i,iterCount=0,domainCount=0,maxIter=50
real*8::a,b,c
real*8::tol=1.d-8

!find interval endpoints
a=-10.
b=10.

print*,'a,b,f(a),f(b)=',a,b,f(a),f(b)

do while(sign(1._8,f(a))==sign(1._8,f(b)) .and. domainCount<20)
  a=2.*a
  b=2.*b
  domainCount=domainCount+1
  print*,'Extending interval: a,b,f(a),f(b)=',a,b,f(a),f(b)
end do

!Compute midpoint
c=0.5*(a+b)

print*,'c,f(c)=',c,f(c)

do while(iterCount<=maxIter) 
  print*,'Iteration#=',iterCount
  !print*,'a,b,c=',a,b,c
  !print*,'f(a),f(b),f(c)=',f(a),f(b),f(c)
  if(sign(1._8,f(a)) .ne. sign(1._8,f(c)))then
    b=c
  end if
  if(sign(1._8,f(b)) .ne. sign(1._8,f(c)))then
    a=c
  end if
  c=0.5*(a+b)

  if(abs(f(c))<tol)then
    print*,'Iterations converged suceesfully.'
    print*,'root=',c
    exit
  end if

  if(abs(f(c))>tol .and. iterCount==40 .and. i<5)then
    print*,'Root finder has not converged yet.'
    print*,'Adding 10 more iterations.'
    maxIter=maxIter+10
    i=i+1
    !exit
  end if
  
  if(abs(f(c))>tol .and. iterCount==40 .and. i==5)then
    print*,'Root finder fails to converge.'
    STOP
  end if
  
  iterCount=iterCount+1

  !print*,'a,b,f(c)=',a,b,f(c)
end do


contains

function f(x) result(fx)
real*8,intent(in)::x

real*8::fx,gnp,dw

gnp=4.
dw=13.


fx=(3.-x)*(dw**(4.-x)-1.)
fx=fx/((4.-x)*(dw**(3.-x)-1.))
fx=fx-gnp
end function f

end program bisectionRoot

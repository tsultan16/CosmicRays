program exponentialTest
implicit none

integer::i

do i=700,750,1
  print*,'exp(',i*1.0,')=',exp(i*1._8)
end do

end program exponentialTest

      subroutine kernel(bcarr,dim,bandwidth,tt,ker)
      implicit none
      real, intent(in), dimension(dim,4) :: bcarr
      real :: ker, kers0, bandwidth, tt, pdf
      integer :: dim, i
      ker=0
      kers0=0
      do i=1,dim
        if (abs(tt-bcarr(i,2)) <= bandwidth) then 
            pdf = (1-((tt-bcarr(i,2))/bandwidth)**2)*3/4
                else 
                    pdf = 0
        end if
      ker = ker + pdf/bandwidth*bcarr(i,3)
      kers0 = kers0 + pdf/bandwidth
      end do
      ker=ker/kers0
      if (kers0==0) then 
           ker=0
      end if
      end subroutine kernel
      
      subroutine kernel2(bcarr,dim,bandwidth,tt,ker)
      implicit none
      real, intent(in), dimension(dim,4) :: bcarr
      real :: ker, kers0, bandwidth, tt, pdf
      integer :: dim, i
      ker=0
      kers0=0
      do i=1,dim
        if (abs(tt-bcarr(i,2)) <= bandwidth) then 
            pdf = (1-abs((tt-bcarr(i,2))/bandwidth))
                else 
                    pdf = 0
        end if
      ker = ker + pdf/bandwidth*bcarr(i,3)
      kers0 = kers0 + pdf/bandwidth
      end do
      ker=ker/kers0
      if (kers0==0) then 
           ker=0
      end if
      end subroutine kernel2
      
       subroutine inter(ball,idall,dim, length,t,tau,yinte)
       implicit none
       real, intent(in), dimension(dim,4) :: ball
       real, intent(in), dimension(length,1) :: idall
       real, intent(inout), dimension(length,1) :: yinte
       real :: t, tau, t1, t2, y1,y2, y
       integer :: dim, length,num,i,j

       do i=1,length
         t1=-1
         t2=tau
         y1=0
         y2=0
         num=0
         y=0
         do j=1,dim
           if (ball(j,1)==idall(i,1)) then
              num=num+1
              y=ball(j,3)
              if (ball(j,2) <= t .AND. ball(j,2) >= t1) then
                 t1=ball(j,2)
                 y1=ball(j,3)
              end if
              if (ball(j,2) >= t .AND. ball(j,2) <= t2) then
                 t2=ball(j,2)
                 y2=ball(j,3)
              end if
            end if
         end do
            if (num==1) then 
               yinte(i,1)=y
            else if (t2==tau) then
               yinte(i,1)=y1 
            else if (t1==-1) then
               yinte(i,1)=y2
            else if (t1==t2) then
               yinte(i,1)=y1
            else
               yinte(i,1)=y1+(t-t1)*(y2-y1)/(t2-t1)
            end if
       end do
       end subroutine inter
       
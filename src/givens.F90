subroutine givens ( a, lda, n, tol )
  !      use FLOPS
  use tools, only: dp
  implicit   none
  integer    lda, n
  real(dp)    a( lda, *), tol

  integer    i,j,k,lastr
  integer(8) count0, count1, count2, count3, count4, count5, count6, &
       count7,count8,count9,count10,count11
  real(dp)    r,x,y,z, b,c,d,e,f,g, t0,t1,t2,t3


  count1 =0 
  count2 =0 
  count3 =0 
  count4 =0 
  count5 =0 
  count6 =0 
  count7 =0 
  count8 =0 
  count9 =0 
  count10 =0 
  count11 =0 
  count0 =0 

  !     zero out upper triangle

  do   j  =  1, n - 1

     !     triple rotations end at multiple of 3

     lastr  =  j + 1 + mod( n - j, 3 )
     do   i  =  n, lastr, - 3 

        count0 = count0 + 1
        ! count0*4
        b =  a( j, i - 1 )
        c =  a( j, i )
        x = b*b + c*c
        if  (  x .gt. tol  )  then       ! test rotation #1

           count1 = count1 + 1
           ! count1*(15+n*18)
           r = x**( -0.5_dp )
           b = b*r
           c = c*r 

           d =  a( j, i - 2 )
           y = d*d + x
           r = y**( -0.5_dp )
           d = d*r
           e = sqrt(x)*r

           f =  a( j, i - 3 )
           z = f*f + y
           r = z**( -0.5_dp )
           f = f*r
           g = sqrt(y)*r

           !     triple-rotation

           !DIR$ ASSUME_ALIGNED A: 64
           !DIR$ ASSUME (mod(lda,8) .eq. 0)
           !DIR$ unroll(8)
           do k = 1,n
              t3 = a( k, i - 3 )
              t2 = a( k, i - 2 )
              t1 = a( k, i - 1 )
              t0 = a( k, i     )
              a( k, i - 3 ) =  t3*f + g*( t2*d + e*( t1*b + t0*c ) )
              a( k, i - 2 ) = -t3*g + f*( t2*d + e*( t1*b + t0*c ) )
              a( k, i - 1 ) = -t2*e + d*( t1*b + t0*c )
              a( k, i     ) = -t1*c + t0*b
           end  do
        else          !  skip rotation 1

           ! (count0-count1)*(4)
           b =  a( j, i - 2 )
           c =  a( j, i - 1 )
           x = b*b + c*c
           if  (  x .gt. tol  )  then       ! test rotation #2
              count3 = count3+1
              ! count3*(9+n*12)
              r = x**( -0.5_dp )
              b = b*r
              c = c*r 

              d =  a( j, i - 3 )
              y = d*d + x
              r = y**( -0.5_dp )
              d = d*r
              e = sqrt(x)*r

              !     double-rotation 2,3
              !DIR$ ASSUME_ALIGNED A: 64
              !DIR$ ASSUME (mod(lda,8) .eq. 0)
              !DIR$ unroll(8)
              do k = 1,n
                 t2 = a( k, i - 3 )
                 t1 = a( k, i - 2 )
                 t0 = a( k, i - 1 )
                 a( k, i - 3 ) =  t2*d + e*( t1*b + t0*c )
                 a( k, i - 2 ) = -t2*e + d*( t1*b + t0*c )
                 a( k, i - 1 ) = -t1*c + t0*b
              end  do
           else            !  skip rotation 2

              !     rotation 3

              ! ((count0-count1) -(count3))*4
              b =  a( j, i - 3 )
              c =  a( j, i - 2 )
              x = b*b + c*c
              if  (  x .gt. tol  )  then       ! test rotation #3
                 count5 = count5 +1
                 ! count5*(3+n*6)
                 r = x**( -0.5_dp )
                 b = b*r
                 c = c*r
                 !DIR$ ASSUME_ALIGNED A: 64
                 !DIR$ ASSUME (mod(lda,8) .eq. 0)
                 !DIR$ unroll(8)
                 do k = 1,n
                    t1 = a( k, i - 3 )
                    t2 = a( k, i - 2 )
                    a( k, i - 3 ) =  t1*b + t2*c
                    a( k, i - 2 ) = -t1*c + t2*b
                 end  do
              end  if         !  test rotation #3
           end  if         !  test rotation #2
        end  if         !  test rotation #1
     end  do         !  columns, i


     !     1 or 2 remainder rotations

     if  (  lastr - j - 1  .eq. 2  )  then
        count6 = count6 + 1
        ! count6*4
        i =  lastr - 1
        b =  a( j, i - 1 )
        c =  a( j, i     )
        x = b*b + c*c
        if  (  x .gt. tol  )  then       ! test rotation #1
           count7 = count7 + 1
           ! count7*(9+n*12)
           r = x**( -0.5_dp )
           b = b*r
           c = c*r 

           d =  a( j, i - 2 )
           y = d*d + x
           r = y**( -0.5_dp )
           d = d*r
           e = sqrt(x)*r

           !     double-rotation 1,2
           !DIR$ ASSUME_ALIGNED A: 64
           !DIR$ ASSUME (mod(lda,8) .eq. 0)
           !DIR$ unroll(8)
           do k = 1,n
              t2 = a( k, i - 2 )
              t1 = a( k, i - 1 )
              t0 = a( k, i     )
              a( k, i - 2 ) =  t2*d + e*( t1*b + t0*c )
              a( k, i - 1 ) = -t2*e + d*( t1*b + t0*c )
              a( k, i     ) = -t1*c + t0*b
           end  do
        else            !  skip rotation 1

           !     rotation 2

           ! (count6-count7)*4
           b =  a( j, i - 2 )
           c =  a( j, i - 1 )
           x = b*b + c*c
           if  (  x .gt. tol  )  then       ! test rotation #2
              count9 = count9 + 1
              ! count9*(3+n*6)
              r = x**( -0.5_dp )
              b = b*r
              c = c*r
              !DIR$ ASSUME_ALIGNED A: 64
              !DIR$ ASSUME (mod(lda,8) .eq. 0)
              !DIR$ unroll(8)
              do k = 1,n
                 t1 = a( k, i - 2 )
                 t2 = a( k, i - 1 )
                 a( k, i - 2 ) =  t1*b + t2*c
                 a( k, i - 1 ) = -t1*c + t2*b
              end  do
           end  if         !  test rotation #2
        end  if         !  test rotation #1

     else  if  (  lastr - j - 1  .eq. 1  )  then

        !     rotation 1
        count10 = count10 + 1
        ! count10*(4)
        i =  lastr - 1
        b =  a( j, i - 1 )
        c =  a( j, i     )
        x = b*b + c*c
        if  (  x .gt. tol  )  then       ! test rotation #1
           count11 = count11 + 1
           !count11*(3+6*n)
           r = x**( -0.5_dp )
           b = b*r
           c = c*r
           !DIR$ ASSUME_ALIGNED A: 64
           !DIR$ ASSUME (mod(lda,8) .eq. 0)
           !DIR$ unroll(8)
           do k = 1,n
              t1 = a( k, i - 1 )
              t2 = a( k, i     )
              a( k, i - 1 ) =  t1*b + t2*c
              a( k, i     ) = -t1*c + t2*b
           end  do
        end  if         !  test rotation #1
     end  if         !  1 or 2 rotations left
  end  do         !  rows, j

  !FLOP = count0*4 + count1*(15+n*18) + (count0-count1)*(4) + count3*(9+n*12) &
  !     + (count0-count1-count3)*4 + count5*( 3 + 6*n) + count6*4 + count7*(9+n*12) &
  !     + (count6-count7)*4 + count9*(3+n*6) + count10*(4) + count11*(3+6*n)
end subroutine givens


subroutine givens_single ( a, lda, n, tol )
  use tools, only : dp
  implicit   none
  integer    lda, n
  real(dp)    a( lda, *), r,cx,sx, tol, t1,t2  !  t0,u0, u1,u2, t3,u3                                                                                                                             
  integer    i,j,k

  !  zero out upper triangle

  do   j = 1,n-1
     do   i = n,j+1,-1

        cx =  a(j,i-1)
        sx =  a(j,i)
        r = cx*cx + sx*sx
        if  (  abs( r ) .gt. tol  )   then
           r = r**( -0.5_dp )
           cx = cx*r
           sx = sx*r

           !   do the matrix multiply
           do k = 1,n
              t1 = a(k,i-1)
              t2 = a(k,i)
              a(k,i-1) = t1*cx + t2*sx
              a(k,i)  = -t1*sx + t2*cx
           end  do

        end  if
     end  do
  end  do

end subroutine givens_single

subroutine givens_orig( a, lda, nord, tol )
  use tools, only: dp
  implicit   none
  integer    lda, nord, nonzero
  !!  unroll,
  !     real*8     a( lda, *), tol, r,cx,sx,t1,t2
  !
  real(dp)    a( lda, *), tol, r,cx,sx, t1, t2 !t0,u0, u1, u2, t3,u3
  integer    nrem

  integer    i,j,k
  logical    converged


  !!  unroll,
  nrem  =  mod( nord, 4 )


  converged = .false.
  do while ( .not. converged )
     nonzero = 0

     do   i = 2, nord
        do   j = 1, i - 1

           !     scan for non-zero elements

           if ( abs( a( j, i ) ) .gt. tol ) then
              nonzero = nonzero + 1

              !     compute sin,cos: Ajj is the diagonal element in the
              !     row-wise rotation (on the right) for the upper-triangle element

              r = ( a( j, i )**2 + a( j, j )**2 )**( -0.5_dp )
              cx =  a( j, j )*r
              sx = -a( j, i )*r

              !     column-update dominates the cost and can be vectorized

              !!  unroll inner loop by hand to expose vectorization
              !
              !!  original code,
              do k = 1, nord
                 t1 = a( k, j )
                 t2 = a( k, i )
                 a( k, j ) = t1*cx - t2*sx
                 a( k, i ) = t1*sx + t2*cx
              end  do
              !
              !!  unrolled,
              !
              !       do k = 1, nord - nrem, 4
              !       t0 = a( k    , j )
              !       t1 = a( k + 1, j )
              !       t2 = a( k + 2, j )
              !       t3 = a( k + 3, j )
              !       u0 = a( k    , i )
              !       u1 = a( k + 1, i )
              !       u2 = a( k + 2, i )
              !       u3 = a( k + 3, i )
              !       a( k    , j ) = t0*cx - u0*sx
              !       a( k + 1, j ) = t1*cx - u1*sx
              !       a( k + 2, j ) = t2*cx - u2*sx
              !       a( k + 3, j ) = t3*cx - u3*sx
              !       a( k    , i ) = t0*sx + u0*cx
              !       a( k + 1, i ) = t1*sx + u1*cx
              !       a( k + 2, i ) = t2*sx + u2*cx
              !       a( k + 3, i ) = t3*sx + u3*cx
              !       end  do
              ! !
              ! !!  remainder loop,
              ! !
              !       do k = nord - nrem + 1, nord
              !       t0 = a( k, j )
              !       u0 = a( k, i )
              !       a( k, j ) = t0*cx - u0*sx
              !       a( k, i ) = t0*sx + u0*cx
              !       end  do

           end  if     !   zero element test
        end  do
     end  do

     converged = nonzero .eq. 0

  end  do     !   lower-triangle zeroed test

  !     the product of the diagonal elements for the
  !     determinant is done in the calling routine
end subroutine givens_orig

#ifdef GIVENS_BGQ
subroutine givens_bgq ( a, lda, n, tol )
  !use FLOPS
  implicit   none
  integer    lda, n
  real(8)    a( lda, *), tol

  integer    i,j,k,lastr,nrem
  integer(8)  count0, count1, count2, count3, count4, count5, count6, &
       count7,count8,count9,count10,count11
  real(8)    r,x,y,z, b,c,d,e,f,g, t0,t1,t2,t3, h,s,u,v
  vector( real(8) ) vt0, vt1, vt2, vt3, vb, vc, vd, ve, vf, vg, vp, vh, vs, vv, vw, vq
  vector( real(8) ) vt00, vt10, vt20, vt30, vp0, vh0, vs0, vv0, vw0, vq0
  vector( real(8) ) vu, vu0

  count0 =0 
  count1 =0 
  count2 =0 
  count3 =0 
  count4 =0 
  count5 =0 
  count6 =0 
  count7 =0 
  count8 =0 
  count9 =0 
  count10 =0 
  count11 =0 

  nrem  =  mod( n, 8 )

  do   j  =  1, n - 1

     lastr  =  j + 1 + mod( n - j, 3 )
     do   i  =  n, lastr, - 3 

        b =  a( j, i - 1 )
        c =  a( j, i )
        x = b*b + c*c
        count0 = count0 + 1
        ! count0*4
        if  (  x .gt. tol  )  then       ! test rotation #1
           count1 = count1 + 1
           ! count1*(21+n/8*36*4)
           ! count1*(21+n*18)
           r = x**( -0.5d0 )
           b = b*r
           c = c*r 

           d =  a( j, i - 2 )
           y = d*d + x
           r = y**( -0.5d0 )
           d = d*r
           e = sqrt(x)*r

           f =  a( j, i - 3 )
           z = f*f + y
           r = z**( -0.5d0 )
           f = f*r
           g = sqrt(y)*r

           vb = VEC_SPLATS( b )
           vc = VEC_SPLATS( c )
           vd = VEC_SPLATS( d )
           ve = VEC_SPLATS( e )
           vf = VEC_SPLATS( f )
           vg = VEC_SPLATS( g )

           do k = 1, n-nrem, 8
              vt3 = VEC_LD( 0, a( k,i-3 ) )
              vt30 = VEC_LD( 0, a( k+4,i-3 ) )

              vt2 = VEC_LD( 0, a( k,i-2 ) )
              vt20 = VEC_LD( 0, a( k+4,i-2 ) )

              vt1 = VEC_LD( 0, a( k,i-1 ) )
              vt10 = VEC_LD( 0, a( k+4,i-1 ) )

              vt0 = VEC_LD( 0, a( k,i ) )
              vt00 = VEC_LD( 0, a( k+4,i ) )

              !IBM* PREFETCH_FOR_LOAD( a( k+16,i-3 ), a( k+16,i-2 ), a( k+16,i-1 ), a( k+16,i ) )

              vp = VEC_MUL( vt0, vb )
              vp0 = VEC_MUL( vt00, vb )

              vp = VEC_NMSUB( vt1, vc, vp )
              vp0 = VEC_NMSUB( vt10, vc, vp0 )

              call VEC_ST( vp, 0, a( k,i ) )
              call VEC_ST( vp0, 0, a( k+4,i ) )

              vh = VEC_MUL( vt1, vb )
              vh0 = VEC_MUL( vt10, vb )

              vh = VEC_MADD( vt0, vc, vh )
              vh0 = VEC_MADD( vt00, vc, vh0 )

              vs = VEC_MUL( vt2, vd )
              vs0 = VEC_MUL( vt20, vd )

              vs = VEC_MADD( ve, vh, vs )
              vs0 = VEC_MADD( ve, vh0, vs0 )

              vv = VEC_MUL( vt3, vf )
              vv0 = VEC_MUL( vt30, vf )

              vv = VEC_MADD( vg, vs, vv )
              vv0 = VEC_MADD( vg, vs0, vv0 )

              call VEC_ST( vv, 0, a( k,i-3 ) )
              call VEC_ST( vv0, 0, a( k+4,i-3 ) )

              vq = VEC_MUL( vt3, vg )
              vq0 = VEC_MUL( vt30, vg )

              vq = VEC_MSUB( vf, vs, vq )
              vq0 = VEC_MSUB( vf, vs0, vq0 )

              call VEC_ST( vq, 0, a( k,i-2 ) )
              call VEC_ST( vq0, 0, a( k+4,i-2 ) )

              vw = VEC_MUL( vt2, ve )
              vw0 = VEC_MUL( vt20, ve )

              vw = VEC_MSUB( vd, vh, vw )
              vw0 = VEC_MSUB( vd, vh0, vw0 )

              call VEC_ST( vw, 0, a( k,i-1 ) )
              call VEC_ST( vw0, 0, a( k+4,i-1 ) )
           end  do

           do k = n - nrem + 1, n
              t3 = a( k, i - 3 )
              t2 = a( k, i - 2 )
              t1 = a( k, i - 1 )
              t0 = a( k, i     )
              a( k, i - 3 ) =  t3*f + g*( t2*d + e*( t1*b + t0*c ) )
              a( k, i - 2 ) = -t3*g + f*( t2*d + e*( t1*b + t0*c ) )
              a( k, i - 1 ) = -t2*e + d*( t1*b + t0*c )
              a( k, i     ) = -t1*c + t0*b
           end  do



        else          !  skip rotation 1

           b =  a( j, i - 2 )
           c =  a( j, i - 1 )
           x = b*b + c*c

           !(count0-count1)*4
           if  (  x .gt. tol  )  then       ! test rotation #2

              count3 = count3+1
              ! count3*(13+n/8*24*4)
              ! count3*(13+n*12)
              r = x**( -0.5d0 )
              b = b*r
              c = c*r 

              d =  a( j, i - 3 )
              y = d*d + x
              r = y**( -0.5d0 )
              d = d*r
              e = sqrt(x)*r

              vb = VEC_SPLATS( b )
              vc = VEC_SPLATS( c )
              vd = VEC_SPLATS( d )
              ve = VEC_SPLATS( e )

              !     double-rotation 2,3
              do k = 1,n-nrem,8
                 !t2 = a( k, i - 3 )
                 vt2 = VEC_LD( 0, a( k,i-3 ) )
                 vt20 = VEC_LD( 0, a( k+4,i-3 ) )

                 !t1 = a( k, i - 2 )
                 vt1 = VEC_LD( 0, a( k,i-2 ) )
                 vt10 = VEC_LD( 0, a( k+4,i-2 ) )

                 !t0 = a( k, i - 1 )
                 vt0 = VEC_LD( 0, a( k,i - 1 ) )
                 vt00 = VEC_LD( 0, a( k+4,i - 1 ) )

                 !IBM* PREFETCH_FOR_LOAD( a( k+12,i-3 ), a( k+12,i-2 ), a( k+12,i - 1) )

                 !a( k, i - 1 ) = t0 * b - t1 * c
                 !u = t0 * b
                 vu = VEC_MUL( vt0, vb )
                 vu0 = VEC_MUL( vt00, vb )

                 !u = u - t1 * c
                 vu = VEC_NMSUB( vt1, vc, vu )
                 vu0 = VEC_NMSUB( vt10, vc, vu0 )

                 !a( k,i - 1 ) = u
                 call VEC_ST( vu, 0, a( k,i - 1 ) )
                 call VEC_ST( vu0, 0, a( k+4,i - 1 ) )

                 !p = t1 * b + t0 * c
                 !p = t1 * b
                 vp = VEC_MUL( vt1, vb )
                 vp0 = VEC_MUL( vt10, vb )

                 !p = p + t0 * c
                 vp = VEC_MADD( vt0, vc, vp )
                 vp0 = VEC_MADD( vt00, vc, vp0 )

                 !u = t2 * d
                 vu = VEC_MUL( vt2,vd)
                 vu0 = VEC_MUL( vt20,vd)

                 !u = u + e * p
                 vu = VEC_MADD( ve, vp, vu )
                 vu0 = VEC_MADD( ve, vp0, vu0 )

                 !a( k, i-3 ) = u
                 call VEC_ST( vu, 0, a( k,i-3 ) )
                 call VEC_ST( vu0, 0, a( k+4,i-3 ) )

                 !a( k, i - 3 ) =  t2*d + e * p
                 !u = d * p
                 vu = VEC_MUL( vd, vp )
                 vu0 = VEC_MUL( vd, vp0 )

                 !u = u - t2 * e
                 vu = VEC_NMSUB( vt2, ve, vu )
                 vu0 = VEC_NMSUB( vt20, ve, vu0 )

                 !a( k,i-2 ) = u
                 call VEC_ST( vu, 0, a( k,i-2 ) )
                 call VEC_ST( vu0, 0, a( k+4,i-2 ) )

                 !a( k, i - 2 ) = d * p - t2 * e
              end  do

              do k = n-nrem+1,n
                 t2 = a( k, i - 3 )
                 t1 = a( k, i - 2 )
                 t0 = a( k, i - 1 )
                 a( k, i - 3 ) =  t2*d + e*( t1*b + t0*c )
                 a( k, i - 2 ) = -t2*e + d*( t1*b + t0*c )
                 a( k, i - 1 ) = -t1*c + t0*b
              end  do



           else            !  skip rotation 2

              !     rotation 3
              b =  a( j, i - 3 )
              c =  a( j, i - 2 )
              x = b*b + c*c
              ! (count0-count1-count3)*4
              if  (  x .gt. tol  )  then       ! test rotation #3
                 r = x**( -0.5d0 )
                 b = b*r
                 c = c*r
                 count5 = count5 +1
                 ! count5 * ( 5 + n/8*12*4)
                 ! count5 * ( 5 + n*6)
                 vb = VEC_SPLATS( b )
                 vc = VEC_SPLATS( c )

                 do k = 1,n-nrem,8

                    !                  t1 = a( k, i - 3 )
                    vt1 = VEC_LD( 0, a( k, i - 3 ) )
                    vt10 = VEC_LD( 0, a( k+4, i - 3 ) )

                    !                  t2 = a( k, i - 2 )
                    vt2 = VEC_LD( 0, a( k, i - 2 ) )
                    vt20 = VEC_LD( 0, a( k+4, i - 2 ) )

                    !IBM* PREFETCH_FOR_LOAD(a( k+16,i-2 ), a( k+16,i -3 ) )

                    vu = VEC_MUL( vt2, vc )
                    vu0 = VEC_MUL( vt20, vc )
                    !                  u = t2*c
                    vu = VEC_MADD( vt1, vb, vu )
                    vu0 = VEC_MADD( vt10, vb, vu0 )
                    !                  u = u + t1*b

                    call VEC_ST( vu, 0, a( k, i - 3 ) )
                    call VEC_ST( vu0, 0, a( k+4, i - 3 ) )
                    !                  a( k, i - 3 ) = u                 
                    !                  a( k, i - 3 ) =  t1*b + t2*c

                    !                  v = t2*b
                    vv = VEC_MUL( vt2, vb )
                    vv0 = VEC_MUL( vt20, vb )
                    !                  v = v - t1*c
                    vv = VEC_NMSUB( vt1, vc, vv )
                    vv0 = VEC_NMSUB( vt10, vc, vv0 )

                    call VEC_ST( vv, 0, a( k, i - 2 ) )
                    call VEC_ST( vv0, 0, a( k+4, i - 2 ) )
                    !                  a( k, i - 2 ) = v
                    !                  a( k, i - 2 ) = -t1*c + t2*b


                 end  do

                 do k = n-nrem+1,n
                    t1 = a( k, i - 3 )
                    t2 = a( k, i - 2 )
                    a( k, i - 3 ) =  t1*b + t2*c
                    a( k, i - 2 ) = -t1*c + t2*b
                 end  do


              end  if         !  test rotation #3
           end  if         !  test rotation #2
        end  if         !  test rotation #1
     end  do         !  columns, i


     !     1 or 2 remainder rotations

     if  (  lastr - j - 1  .eq. 2  )  then
        i =  lastr - 1
        b =  a( j, i - 1 )
        c =  a( j, i     )
        x = b*b + c*c
        count6 = count6 + 1
        ! count6*(4)
        if  (  x .gt. tol  )  then       ! test rotation #1
           r = x**( -0.5d0 )
           b = b*r
           c = c*r 

           d =  a( j, i - 2 )
           y = d*d + x
           r = y**( -0.5d0 )
           d = d*r
           e = sqrt(x)*r

           count7 = count7 + 1
           ! count7*( 13 + n/8 * 24 * 4)
           ! count7*( 13 + n*12)
           vb = VEC_SPLATS( b )
           vc = VEC_SPLATS( c )
           vd = VEC_SPLATS( d )
           ve = VEC_SPLATS( e )


           !     double-rotation 1,2
           do k = 1,n-nrem,8
              vt2 = VEC_LD( 0, a( k,i-2 ) )
              vt20 = VEC_LD( 0, a( k+4,i-2 ) )

              vt1 = VEC_LD( 0, a( k,i-1 ) )
              vt10 = VEC_LD( 0, a( k+4,i-1 ) )

              vt0 = VEC_LD( 0, a( k,i ) )
              vt00 = VEC_LD( 0, a( k+4,i ) )

              !IBM* PREFETCH_FOR_LOAD( a( k+12,i-2 ), a( k+12,i-1 ), a( k+12,i ) )

              vu = VEC_MUL( vt0, vb )
              vu0 = VEC_MUL( vt00, vb )

              vu = VEC_NMSUB( vt1, vc, vu )
              vu0 = VEC_NMSUB( vt10, vc, vu0 )

              call VEC_ST( vu, 0, a( k,i ) )
              call VEC_ST( vu0, 0, a( k+4,i ) )

              vp = VEC_MUL( vt1, vb )
              vp0 = VEC_MUL( vt10, vb )

              vp = VEC_MADD( vt0, vc, vp )
              vp0 = VEC_MADD( vt00, vc, vp0 )

              vu = VEC_MUL( vt2,vd)
              vu0 = VEC_MUL( vt20,vd)

              vu = VEC_MADD( ve, vp, vu )
              vu0 = VEC_MADD( ve, vp0, vu0 )

              call VEC_ST( vu, 0, a( k,i-2 ) )
              call VEC_ST( vu0, 0, a( k+4,i-2 ) )

              vu = VEC_MUL( vd, vp )
              vu0 = VEC_MUL( vd, vp0 )

              vu = VEC_NMSUB( vt2, ve, vu )
              vu0 = VEC_NMSUB( vt20, ve, vu0 )

              call VEC_ST( vu, 0, a( k,i-1 ) )
              call VEC_ST( vu0, 0, a( k+4,i-1 ) )

           end  do

           do k = n - nrem + 1,n
              t2 = a( k, i - 2 )
              t1 = a( k, i - 1 )
              t0 = a( k, i     )
              a( k, i - 2 ) =  t2*d + e*( t1*b + t0*c )
              a( k, i - 1 ) = -t2*e + d*( t1*b + t0*c )
              a( k, i     ) = -t1*c + t0*b
           end  do



        else            !  skip rotation 1

           !     rotation 2
           b =  a( j, i - 2 )
           c =  a( j, i - 1 )
           x = b*b + c*c
           ! (count6-count7)*(4)
           if  (  x .gt. tol  )  then       ! test rotation #2
              r = x**( -0.5d0 )
              b = b*r
              c = c*r

              count9 = count9 + 1
              !count9*( 5+ n/8 * 12 *4)
              !count9*( 5+ n*6)
              vb = VEC_SPLATS( b )
              vc = VEC_SPLATS( c )

              do k = 1,n-nrem,8

                 !               t1 = a( k, i - 2 )
                 vt1 = VEC_LD( 0, a(k,i-2) )
                 vt10 = VEC_LD( 0, a(k+4,i-2) )

                 !               t2 = a( k, i - 1  )
                 vt2 = VEC_LD( 0, a(k,i - 1 ) )
                 vt20 = VEC_LD( 0, a(k+4,i - 1) )

                 !IBM* PREFETCH_FOR_LOAD(a( k+16,i-2 ), a( k+16,i - 1) )
                 !               u = t2*c
                 vu = VEC_MUL(vt2, vc )
                 vu0 = VEC_MUL(vt20, vc )
                 !               u = u + t1*b
                 vu = VEC_MADD( vt1, vb, vu )
                 vu0 = VEC_MADD( vt10, vb, vu0 )

                 call VEC_ST( vu, 0, a( k, i - 2 ) )
                 call VEC_ST( vu0, 0, a( k+4, i - 2 ) )
                 !               a( k, i - 2 ) =  u
                 !               a( k, i - 2 ) =  t1*b + t2*c

                 !               v = t2*b
                 vp = VEC_MUL( vt2, vb)
                 vp0 = VEC_MUL( vt20, vb)
                 !               v = v - t1*c
                 vp = VEC_NMSUB( vt1, vc, vp )
                 vp0 = VEC_NMSUB( vt10, vc, vp0 )

                 call VEC_ST( vp, 0, a( k, i - 1 ) )
                 call VEC_ST( vp0, 0, a( k+4, i - 1 ) )
                 !               a( k, i - 1  ) = v
                 !               a( k, i - 1 ) = -t1*c + t2*b
              end  do


              do k = n - nrem + 1,n
                 t1 = a( k, i - 2 )
                 t2 = a( k, i - 1 )
                 a( k, i - 2 ) =  t1*b + t2*c
                 a( k, i - 1 ) = -t1*c + t2*b
              end  do


           end  if         !  test rotation #2
        end  if         !  test rotation #1

     else  if  (  lastr - j - 1  .eq. 1  )  then

        !     rotation 1
        count10 = count10 + 1
        ! count10*4
        i =  lastr - 1
        b =  a( j, i - 1 )
        c =  a( j, i     )
        x = b*b + c*c
        if  (  x .gt. tol  )  then       ! test rotation #1
           count11 = count11 + 1
           ! count11*(5 + n/8 * 12 * 4)
           ! count11*(5 + n *6 )
           r = x**( -0.5d0 )
           b = b*r
           c = c*r

           vb = VEC_SPLATS( b )
           vc = VEC_SPLATS( c )

           do k = 1,n-nrem,8

              !               t1 = a( k, i - 1 )
              vt1 = VEC_LD( 0, a(k,i-1) )
              vt10 = VEC_LD( 0, a(k+4,i-1) )

              !               t2 = a( k, i     )
              vt2 = VEC_LD( 0, a(k,i) )
              vt20 = VEC_LD( 0, a(k+4,i) )

              !IBM* PREFETCH_FOR_LOAD(a( k+16,i-1 ), a( k+16,i ) )
              !               u = t2*c
              vu = VEC_MUL(vt2, vc )
              vu0 = VEC_MUL(vt20, vc )
              !               u = u + t1*b
              vu = VEC_MADD( vt1, vb, vu )
              vu0 = VEC_MADD( vt10, vb, vu0 )

              call VEC_ST( vu, 0, a( k, i - 1 ) )
              call VEC_ST( vu0, 0, a( k+4, i - 1 ) )
              !               a( k, i - 1 ) =  u
              !               a( k, i - 1 ) =  t1*b + t2*c

              !               v = t2*b
              vp = VEC_MUL( vt2, vb)
              vp0 = VEC_MUL( vt20, vb)
              !               v = v - t1*c
              vp = VEC_NMSUB( vt1, vc, vp )
              vp0 = VEC_NMSUB( vt10, vc, vp0 )

              call VEC_ST( vp, 0, a( k, i     ) )
              call VEC_ST( vp0, 0, a( k+4, i     ) )
              !               a( k, i     ) = v
              !               a( k, i     ) = -t1*c + t2*b
           end  do

           do k = n-nrem+1,n
              t1 = a( k, i - 1 )
              t2 = a( k, i     )
              a( k, i - 1 ) =  t1*b + t2*c
              a( k, i     ) = -t1*c + t2*b
           end  do


        end  if         !  test rotation #1
     end  if         !  1 or 2 rotations left
  end  do         !  rows, j

  !FLOP = count0*4 + count1*(21+36*4*n/8)
  ! FLOP = FLOP+ (count0*4 + count1*(21+n*18) + (count0-count1)*(4) + count3*(13+n*12) &
  !      + (count0-count1-count3)*(4) + count5*(5 + 6*n) + count6*4 + count7*(13+n*12) &
  !      + (count6-count7)*4 + count9*(5+n*6) + count10*(4) + count11*(5+6*n))

  !print *, "count3", count0, count1, count0-count1, count3, (count0-count1-count3), count5, count6, count7, count6-count7, count9, count10, count11

end subroutine givens_bgq
#endif

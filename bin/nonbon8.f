c     main energy-calc function for OLD pairlists, grid uses nonbond.h !!!
      subroutine nonbon8(iab,xl,xr,fl,fr,wel,wer,chair,chail,ac,rc,
     1 emin,rmin2,iacir,iacil,nonr,nonl,ipon,nonp,
     2 potshape, cdie, swi_on,swi_off, enon,epote)
      implicit none

c     Parameters
      include "max.fin"
      integer iab,nonp,potshape
      integer cdie
      real swi_on, swi_off
      real*8 xl,xr,fl,fr,wel,wer,chair,chail,ac,rc
      real*8 emin,rmin2,enon,epote      
      integer iacir,iacil,ipon,nonl,nonr
      dimension nonr(maxmolpair),nonl(maxmolpair),ipon(99,99),
     1  iacir(maxatom), iacil(maxatom)
      dimension chair(maxatom),chail(maxatom),wer(maxatom),
     1 wel(maxatom),fl(maxatom),fr(maxatom),ac(99,99),rc(99,99),
     2 rmin2(99,99)
      dimension xl(maxatom),xr(maxatom),emin(99,99)

c     Local variables
      real*8 dx,xnull,r2,rr1,rr2,rrd,rr23,rep,vlj,et,charge,alen,
     1 rlen,fb,fdb, rr2a
      real*8 fswi, r, shapedelta
      integer k,ik,i,j,ii,jj,it,jt,ivor
      dimension dx(3)
      real*8 e_min
      
      xnull=0.0d0
      enon=xnull
      epote=xnull
      r2=xnull
      Do 100 ik = 1, nonp ! for all nonbonded-pairs
         i = nonr (ik)
         j = nonl (ik)
         it = iacir (i)! i-atomtype = interaction-atom_type_receptor(i)
         jt = iacil (j)! j-...= ..._ligand
         ii = 3 * (i-1)! i for coords?
         jj = 3 * (j-1)
      !Attractive LJ-ENergy ? or attraction-length? -> in any case: zaehler im LJ-Bruch, über radius
         alen = wel (j) * wer (i) * ac (it, jt)! = ?*? * attractive-X(Atomtype-pair) > aus 1. 99*99 inmput-matrix
      ! Repulsive ...
         rlen = wel (j) * wer (i) * rc (it, jt) ! repulive matrix?
      ! minimum of LJ-Pot.
         e_min = wel (j) * wer (i) * emin (it, jt) ! depth-matrix
      ! switch,makes potential non-LJ but repulsive only (atomtypes just dont go well together)
         ivor = ipon (it, jt)
         charge = wel (j) * wer (i) * chair (i) * chail (j)! chair= charge receptor, chail = ..lig
         r2 = xnull ! useless? already above...
      ! calc square of distance -> add coord-difference squares...
         Do 120 k = 1, 3
            dx (k) = xl (jj+k) - xr (ii+k)
            r2 = r2 + dx (k) ** 2
120      Continue
         If (r2 .Lt. 0.001d0) r2 = 0.001d0 ! if too close to 0: set slightly above -> no 0-division..

         fswi = 1.0d0
         If (swi_on .Gt. 0 .Or. swi_off .Gt. 0) Then
            If (r2 .Ge. (swi_on*swi_on)) Then
               If (r2 .Ge. (swi_off*swi_off)) Then
                  fswi = 0.0d0
               Else
                  r = Sqrt (r2)
            fswi = 1.0d0 - (r-swi_on) / (swi_off-swi_on)
               End If
            End If
         End If
      ! inverse of radius^2
         rr2 = 1.0d0 / r2
         Do 125 k = 1, 3
            dx (k) = rr2 * dx (k)
125      Continue
         et = xnull
         If (charge .Gt. 0.001 .Or. charge .Lt.-0.001) Then
            If (cdie .Eq. 1) Then
               rr1 = 1.0d0 / Sqrt (r2) - 1.0 / 50.0
c            (cap all distances at 50 A)
               If (rr1 .Lt. 0) Then
                  rr1 = 0
               End If
               et = charge * rr1
c       write(*,*), sqrt(r2), et, epote
            Else
               rr2a = rr2 - (1.0/50.0) * (1.0/50.0)
               If (rr2a .Lt. 0) Then
                  rr2a = 0
               End If
c      (cap all distances at 50 A)
               et = charge * rr2a
            End If
            epote = epote + fswi * et
            If (iab .Eq. 1) Then
               Do 130 k = 1, 3
                  If (cdie .Eq. 1) Then
                     If (rr1 .Le. 0) Then
                        fdb = fswi * et * dx (k)
                     Else
                        fdb = fswi * charge * (rr1+1.0/50.0) * dx (k)
                     End If
                  Else
                     If (rr2a .Le. 0) Then
                        fdb = fswi * 2.0d0 * et * dx (k)
                     Else
                        fdb = fswi * 2.0d0 * charge * rr2 * dx (k)
                     End If
                  End If
                  ! force-lig, force-rec
                  fl (jj+k) = fl (jj+k) + fdb
                  fr (ii+k) = fr (ii+k) - fdb
130            Continue
            End If
         End If
         ! 1/r⁶
         rr23 = rr2 ** 3
         If (potshape .Eq. 8) Then
            rrd = rr2
            shapedelta = 2.0D0
         Else If (potshape .Eq. 12) Then
            rrd = rr23
            shapedelta = 6.0D0
         End If

         rep = rlen * rrd
         vlj = (rep-alen) * rr23 !LJ-pot finally calculated
         If (r2 .Lt. rmin2(it, jt)) Then ! is distance under minimum?
            enon = enon + fswi * (vlj+(ivor-1)*e_min) ! if ivor=1 (normal potential) -> last term is 0, else (rep0ulsive-only potential) it shift the pot by 2*e_min
c       write(*,*)'pair',i,j,it,jt,r2,vlj+(ivor-1)*e_min,e_min
c      write(*,*)'pair',i,j,it,jt,r2,vlj+(ivor-1)*emin(it,jt),et,
c     1 emin(it,jt)
            If (iab .Eq. 1) Then
               fb = 6.0D0 * vlj + shapedelta * (rep*rr23)
               Do 135 k = 1, 3
                  fdb = fswi * fb * dx (k)
                  fl (jj+k) = fl (jj+k) + fdb
                  fr (ii+k) = fr (ii+k) - fdb
135            Continue
            End If
         Else
            enon = enon + fswi * ivor * vlj 
c      write(*,*)'pair',i,j,it,jt,r2,ivor*vlj,et,
c     1 emin(it,jt)
            If (iab .Eq. 1) Then
               fb = 6.0D0 * vlj + shapedelta * (rep*rr23)

               Do 145 k = 1, 3
                  fdb = fswi * ivor * fb * dx (k)
                  fl (jj+k) = fl (jj+k) + fdb
                  fr (ii+k) = fr (ii+k) - fdb
145            Continue
            End If
         End If
100   Continue
c      write(ERROR_UNIT, *) "nonbon8 ", enon, epote

      Return
      End

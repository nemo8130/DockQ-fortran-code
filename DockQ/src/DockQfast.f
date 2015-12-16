*************************************************************************
!========================================================================
!
!     Calculate LRMS, IRMS, FNAT based on sequence dependent
!     superposition
!     Use all (N, CA, C) backbone atoms 
!     Usage: ./program MODEL.pdb TARGET.pdb T1onT2_irms.pdb T1onT2_lrms.pdb
!     
!========================================================================


      program calRMSD
      PARAMETER(nmax=100000)
      
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)

      character*500 fnam,pdbnat,pdbmod,outname1,outname2,cho,cht,pdb(2)
      character*3 aa(-1:20),seqA(nmax),seqB(nmax)
      character(3)::atomA(nmax),atomB(nmax),res1,res2,res1m,res2m
      character*500 s,du
      character(20)::capstat
      character seq1A(nmax),seq1B(nmax),ali(nmax)
      character sequenceA(nmax),sequenceB(nmax),sequenceM(nmax)
      character(1)::chA(nmax),chB(nmax),chain1,chain2,chainsup,
     &chaincal,chAB(nmax),chstA(nmax),chstB(nmax)

      dimension L_ini(100),iq(nmax)
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10
      dimension xa(nmax),ya(nmax),za(nmax)
      real::lrms,lrms1,lrms2,irms
      integer::nintA(nmax),nintB(nmax),nintAfnat(nmax),nintBfnat(nmax)
      integer::nintAB(nmax)
      real::xat(30),yat(30),zat(30),xbt(30),ybt(30),zbt(30)
      real::xatm(30),yatm(30),zatm(30),xbtm(30),ybtm(30),zbtm(30)
      real::DockQ,c1,c2,c3

ccc   RESIDUE_POINTERS
      integer::istore1(nmax),istore2(nmax),iresA(nmax),iresB(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

      cutdist = 5.00 ! (Between any two atoms)
      cutdistI = 10.00 ! (Between any two atoms)
      cutarb = 16.00 ! (Upper threshold beyond which exhaustive pairwise atomic contact calculation is not considered for a given residue pair)


*****instructions ----------------->
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
        write(*,*)'-------------------------------------------'
        write(*,*)
        write(*,*)'Usage: ./DockQ.exe model.pdb native.pdb'
        write(*,*)'PDB files must be pre-aligned'
        write(*,*)
        write(*,*)'-------------------------------------------'
        goto 9999
      endif
******* options ----------->
      narg=iargc()
      if(narg<2)then
         write(*,*)'Error, please enter two PDB files'
         goto 9999
      endif
      i=0
      j=0
      m_out=0
 115  continue
      i=i+1
      call getarg(i,fnam)
         j=j+1
         pdb(j)=fnam
      if(i.lt.narg)goto 115

      pdbmod = pdb(1)
      pdbnat = pdb(2)

!      call getarg(3,outname1)
!      call getarg(4,outname2)
!      print*,outname1,outname2
******* options ----------->
ccccccccc read data from first CA file:
!==============================================================
!  READ ALL ATOMS FROM pdb file 1 (MODEL)
!==============================================================
      open(unit=10,file=pdbmod,status='old')
      i=0
      ind1=0
 101  read(10,104,end=102) s
      if(s(1:3).eq.'TER') goto 112
      if(s(1:3).eq.'END') goto 102
      if(s(1:4).eq.'ATOM')then
!         if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
!     &        eq.'  CA')then
         if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
            if (s(14:14).ne.'H'.and.s(27:27).ne.'X')then
            i=i+1
            read(s,103)du,atomA(i),seqA(i),chA(i),nresA(i),du,xa(i),
     &ya(i),za(i)
!            write(*,103)du,atomA(i),seqA(i),chA(i),nresA(i),du,xa(i),
!     &ya(i),za(i)
                        if (nresA(i) /= nresA(i-1))then
                        ind1 = ind1 + 1
                        istore1(ind1) = i
                        iresA(ind1) = nresA(i)
                        chstA(ind1) = chA(i)
                        endif
                        
         endif
         endif
      endif
112   continue
      goto 101
 102  continue
 103  format(A13,A3,1x,A3,1X,A1,i4,A4,3F8.3)
 104  format(A100)
      close(10)
      nseqA=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     
!      print*,'chains read from pdb1:  ',chA(1),'  ',chA(nseqA)

ccccccccc read data from second CA file:
!==============================================================
!  READ ALL ATOMS FROM pdb file 2 (TARGET)
!==============================================================
      open(unit=10,file=pdbnat,status='old')
      i=0
      ind2=0
 201  read(10,204,end=202) s
      if(s(1:3).eq.'TER') goto 212
      if(s(1:3).eq.'END') goto 202
      if(s(1:4).eq.'ATOM')then
!         if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
!     &        eq.'  CA')then
         if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
            if (s(14:14).ne.'H'.and.s(27:27).ne.'X')then
            i=i+1
            read(s,203)du,atomB(i),seqB(i),chB(i),nresB(i),du,xb(i),
     &yb(i),zb(i)
!            write(*,203)du,atomB(i),seqB(i),chB(i),nresB(i),du,xb(i),
!     &yb(i),zb(i)
                        if (nresB(i) /= nresB(i-1))then
                        ind2 = ind2 + 1
                        istore2(ind2) = i
                        iresB(ind2) = nresB(i)
                        chstB(ind2) = chB(i)
!            write(*,203)du,atomB(i),seqB(i),chB(i),nresB(i),du,xb(i),
!     &yb(i),zb(i)
                        endif
         endif
         endif
      endif
 212  continue
      goto 201
 202  continue
 203  format(A13,A3,1x,A3,1X,A1,i4,A4,3F8.3)
 204  format(A100)
      close(10)
      nseqB=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
!      print*,'chains read from pdb2:  ',chB(1),'  ',chB(nseqB)
!      print*,'Number of atoms read:',nseqA,'(MODEL) and ',nseqB,
!     &'(NATIVE)'

      if ((chA(1).eq.chB(1)).and.(chA(nseqA).eq.chB(nseqB)))then
      chain1 = chA(1)
      chain2 = chA(nseqA)
      else 
!      print*, chA(1),'  ',chB(1),'      ',chA(nseqA),'  ',chB(nseqB)
      print*,'CHIAN IDs DOES NOT MATCH IN THE TWO FILES. EXITING'
      goto 555
      endif

!      print*,'chain1:', chain1
!      print*,'chain2:', chain2

!      print*,ind1,ind2

!================================ FOR THE LAST RESIDUE ==============================

      istore1(ind1+1) = nseqA
      istore2(ind2+1) = nseqB

!================================================================================
!     DETECT LIGAND (SMALLER PARTNER) AND RECEPTOR (LARGER PARTNER) IN THE NATIVE
!     ASSUME IT WOULD BE THE SAME IN THE MODEL
!     IN CASE OF IDENTICAL SIZE, CONSIDER THE FIRST CHAIN AS RECEPTOR
!================================================================================

      ich1 = 0
      ich2 = 0

      do i = 1,nseqA
           if (chA(i).eq.chain1)then
           ich1 = ich1 + 1    
           elseif (chA(i).eq.chain2)then
           ich2 = ich2 + 1
           endif
      enddo

!      print*,chain1,': ',ich1,'   ',chain2,': ',ich2
!      print*,nseqA,nseqB

      if (ich1 >= ich2)then
!      print*,chain2,' is the Ligand and ',chain1,' is the receptor'
      chainsup = chain1
      chaincal = chain2
      else
!      print*,chain1,' is the Ligand and ',chain2,' is the receptor'
      chainsup = chain2
      chaincal = chain1
      endif

! chainsup:  receptor, chaincal: ligand

      print*,'*********************************************************'
      print*,'*                       DockQ                           *'
      print*,'*   Scoring function for protein-protein docking models *'
      print*,'*   Statistics on CAPRI data:                           *'
      print*,'*    0    <  DockQ <  0.23 - Incorrect                  *'
      print*,'*    0.23 <= DockQ <  0.49 - Acceptable quality         *'
      print*,'*    0.40 <= DockQ <  0.79 - Medium quality             *'
      print*,'*            DockQ >= 0.79 - High quality               *'
      print*,'*   Reference: Sankar Basu and Bjorn Wallner, DockQ:... *'
      print*,'*   For comments, please email: bjornw@ifm.liu.se       *'
      print*,'*********************************************************'

!===================================================================================================
!===================================================================================================
! NATIVE
!===================================================================================================
!===================================================================================================

!      write(27,829)'load pdb ',pdbnat
!      write(27,*)'wireframe off'
!      write(27,*)'ribbon'
!      write(27,*)'color chain'


!      print*,ind1,ind2

!      print*,'Native interfacial conatcts:'

      icontN=0
      iint = 0
      itot=0

      do i = 1,ind2             ! on the index array for the first time
      iresflag1 = 0   ! residue flag
           do j = (i+1),ind2        ! on the index array for the second time
           i1 = 0
           j1 = 0
           iflag1=0
             if ((chstB(i).eq.chain1) .and. (chstB(j).eq.chain2))then  ! receptor vs. ligand contacts
                do ii = istore2(i),istore2(i+1)-1
                i1 = i1 + 1
                xat(i1) = xb(ii)
                yat(i1) = yb(ii)
                zat(i1) = zb(ii)
                res1=seqB(ii)
                enddo
                do jj = istore2(j),istore2(j+1)-1
                j1 = j1 + 1
                xbt(j1) = xb(jj)
                ybt(j1) = yb(jj)
                zbt(j1) = zb(jj)
                res2=seqB(jj)
                enddo
!              print*,iresB(i),chstB(i),i1,iresB(j),chstB(j),j1
!===================================================================================================
!===================================================================================================
!===================================================================================================
!============= CALCULATE CONTACTS BETWEEN iresB(i)-chstB(i) and iresB(j)-chstB(j) ==================
!===================================================================================================
!===================================================================================================
!===================================================================================================
                      do k1 = 1,i1
                          do k2 = 1,j1
                dist=10000.0
                call caldist(xat(k1),yat(k1),zat(k1),xbt(k2),ybt(k2),
     &zbt(k2),dist)
!================================================================= 
! If any pair of atom is more than an arbotary threshold, no point in doing an exhaustive search for that particular pair of residue
!=================================================================
!                               if (dist > cutarb)then
!                               goto 1661
!                               endif
!==============================================================================
! Hit to just the first atomic contact is sufficient to count a res-res contact
!==============================================================================
                               if (dist <= cutdistI)then      ! for IRMS
                               iresflag1=1
                                   if (dist <= cutdist)then   ! for Fnat
                                   iflag1=1
                                   goto 1661
                                   endif
                               endif
                          enddo
                      enddo
1661   continue
!===================================================================================================
!===================================================================================================
!===================================================================================================
             endif

             if (iflag1==1)then
             icontN = icontN + 1                ! iresB(i)-chstB(i) and iresB(j)-chstB(j) are in contact
             nintAfnat(icontN) = iresB(i)
             nintBfnat(icontN) = iresB(j)

!             print*,iresA(i),chstA(i),i1,iresA(j),chstA(j),j1,icontM
!             print*,nintAfnat(icontN),chstA(i),
!     &nintBfnat(icontN),chstB(j)
             endif
           itot = itot+1
           enddo
           if (iresflag1==1)then
           iint = iint + 1
           nintAB(iint) = iresB(i)
           chAB(iint) = chstB(i)
           endif
      enddo

!      print*,'Number of native contacts at the receptor-ligand interface
!     & :',icontN

!      print*,'Number of interfacial residues from the receptor 
!     &defined at 10 A:',iint

!      goto 9999

!============ detect ligand interface residues ===========================

      do i = 1,ind2             ! on the index array for the first time
      iresflag1 = 0   ! residue flag
           do j = 1,ind2        ! on the index array for the second time
           i1 = 0
           j1 = 0
           iflag1=0
!       print*,chstB(i)
            if ((chstB(i).eq.chain2) .and. (chstB(j).eq.chain1))then
! ligand vs. receptor contacts (revert the order)
                do ii = istore2(i),istore2(i+1)-1
                i1 = i1 + 1
                xat(i1) = xb(ii)
                yat(i1) = yb(ii)
                zat(i1) = zb(ii)
                res1=seqB(ii)
                enddo
                do jj = istore2(j),istore2(j+1)-1
                j1 = j1 + 1
                xbt(j1) = xb(jj)
                ybt(j1) = yb(jj)
                zbt(j1) = zb(jj)
                res2=seqB(jj)
                enddo
!===================================================================================================
!===================================================================================================
                      do k1 = 1,i1
                          do k2 = 1,j1
                call caldist(xat(k1),yat(k1),zat(k1),xbt(k2),ybt(k2),
     &zbt(k2),dist)
!================================================================= 
!                               if (dist > cutarb)then
!                               goto 3661
!                               endif
!==============================================================================
                               if (dist <= cutdistI)then      ! for IRMS
                               iresflag1=1
                               goto 3661
                               endif
                          enddo
                      enddo
3661   continue
!===================================================================================================
!===================================================================================================
             endif
           enddo
           if (iresflag1==1)then
           iint = iint + 1
           nintAB(iint) = iresB(i)
           chAB(iint) = chstB(i)
           endif
        enddo


!===================================================================
!     PULL UP ALL NATIVE INTERFACIAL RESIDUES 
!===================================================================

!      print*,'Nintres receptor, ligand defined at 10 A:',iint

!      do i = 1,iint
!      write(*,*)nintAB(i),'  ',chAB(i)
!      enddo

!      print*,'Total number of int res',iint

!===================================================================================================
!===================================================================================================
!  MODEL 
!===================================================================================================
!===================================================================================================

      icontM=0
      ipres=0

      do i = 1,ind1             ! on the index array for the first time
           do j = (i+1),ind1        ! on the index array for the second time
           i1 = 0
           j1 = 0
             if ((chstA(i).eq.chain1) .and. (chstA(j).eq.chain2))then   ! receptor vs. ligand contacts
                do ii = istore1(i),istore1(i+1)-1
                i1 = i1 + 1
                xat(i1) = xa(ii)
                yat(i1) = ya(ii)
                zat(i1) = za(ii)
                enddo
                do jj = istore1(j),istore1(j+1)-1
                j1 = j1 + 1
                xbt(j1) = xa(jj)
                ybt(j1) = ya(jj)
                zbt(j1) = za(jj)
                enddo
!              print*,iresB(i),chstB(i),i1,iresB(j),chstB(j),j1
!===================================================================================================
!===================================================================================================
!===================================================================================================
!============= CALCULATE CONTACTS BETWEEN iresB(i)-chstB(i) and iresB(j)-chstB(j) ==================
!===================================================================================================
!===================================================================================================
!===================================================================================================
                 iflag1=0
                      do k1 = 1,i1
                          do k2 = 1,j1
                call caldist(xat(k1),yat(k1),zat(k1),xbt(k2),ybt(k2),
     &zbt(k2),dist)
!================================================================= 
! If any pair of atom is more than an arbotary threshold, no point in
! doing an exhaustive search for that particular pair of residue
!=================================================================
!                               if (dist > cutarb)then
!                               goto 2661
!                               endif
!==============================================================================
! Hit to just the first atomic contact is sufficient to count a res-res
! contact
!==============================================================================
                               if (dist <= cutdist)then
                               icontM = icontM + 1
                               iflag1=1
                               goto 2661
                               endif
                          enddo
                      enddo
2661   continue
!===================================================================================================
!===================================================================================================
!===================================================================================================
             endif
             if (iflag1==1)then
!             print*,iresA(i),chstA(i),i1,iresA(j),chstA(j),j1,icontM
!             print*,nintAfnat(icontM),chstA(i),
!     &nintBfnat(icontM),chstB(j)
!============================== FIND MATCH TO NATIVE ===============================================              
                do ij = 1,icontN
                  if ((nintAfnat(ij)==iresA(i).and.chstA(i).eq.chain1)
     &.and.((nintBfnat(ij)==iresA(j).and.chstA(j).eq.chain2)))then
!                   print*,nintAfnat(ij),chain1,nintBfnat(ij),chain2
                   ipres=ipres+1
                   endif
                enddo

             endif
           enddo
      enddo

!      print*,icontM,ipres

!      print*,'Native:',icontN,'Model:',icontM,'Preserved:',ipres

         if (icontN == 0)then
         print*,'No Interfacial Contacts found in the native'
         write(245,*)'No Interfacial Contacts found in the native'
         print*,'Exiting'
         goto 9999
         endif


      fnat=float(ipres)/float(icontN)
      inonat=icontM-ipres
      fnonat=float(inonat)/float(icontM)

      write(*,181)'fnat: ',fnat,ipres,'correct out of',icontN,'native 
     &contacts'

      write(*,182)'fnonat: ',fnonat,inonat,'non-native out of',icontM,
     &'model contacts'

181    format(a10,f8.3,':',i5,a20,i5,a20)
182    format(a10,f8.3,':',i5,a20,i5,a20)

!      print*,'------------------------------'


!      goto 999

!=====================================================================================


829   format(a9,a80)

99    format(a7,i3,a6,a1,a5,'*',a3)
856   format(a9)
91    format(i4,2x,a1,5x,i4,2x,a1,2x,i5)
92    format(i3,a1,a3,a1,a1,2x,i3,a1,a3,a1,a1)
93    format(a3,2x,i4,2x,a1,5x,i4,2x,a1,2x,i5)
95    format(i3,a1,a3,a1,a1)

!      print*,'Total number of int res',iint

******************************************************************
*     pickup the common residues:
******************************************************************
      k=0
      do i=1,nseqA
         do j=1,nseqB
            if((nresA(i).eq.nresB(j)) .and. (chA(i).eq.chB(j)).and. 
     &(seqA(i).eq.seqB(j)) .and. (atomA(i).eq.atomB(j)))then
!==================================================================================
!           BACKBONE ATOMS ONLY 
!==================================================================================
             if ((atomA(i).eq.'N  '.or.atomA(i).eq.'CA '
     &.or.atomA(i).eq.'C  '.or.atomA(i).eq.'O  ').and.
     &(atomB(j).eq.'N  '.or.atomB(j).eq.'CA '.or.
     &atomB(j).eq.'C  '.or.atomB(j).eq.'O  '))then
               k=k+1
               iA(k)=i
               iB(k)=j
!       write(191,116)nresA(i),seqA(i),chA(i),atomA(i),'<=>',
!     &nresB(j),seqB(j),chB(j),atomB(j)
!               goto 205
             endif
            endif
         enddo
! 205     continue
      enddo
      n_ali=k                   !number of aligned atoms
      if(n_ali.lt.1)then
        write(*,*)'There is no common residues in the input structures'
        goto 9999
      endif

!      print*,'Number of Aligned Backbone Atoms: ',n_ali

116   format(i4,2x,a3,2x,a1,2x,a4,2x,a3,i4,2x,a3,2x,a1,2x,a4)

***   IRMSD calculation ---------------------------------->

!      print*,'IRMSD:'

      do i = 1,3
      t(i) = 0.000
           do j = 1,3
           u(i,j) = 0.000
           enddo
      enddo

!==================================================================================
      in1 = 0

      do m=1,n_ali
           do i = 1,iint
!==================================================================================
!          INTERFACE ATOMS ONLY
!==================================================================================
           if(nresA(iA(m))==nintAB(i).and.chA(iA(m)).eq.chAB(i))then
           in1 = in1 + 1
           r_1(1,in1)=xa(iA(m))
           r_1(2,in1)=ya(iA(m))
           r_1(3,in1)=za(iA(m))
           r_2(1,in1)=xb(iB(m))
           r_2(2,in1)=yb(iB(m))
           r_2(3,in1)=zb(iB(m))
!         write(*,718)nresA(iA(m)),seqA(iA(m)),atomA(iA(m)),chA(iA(m)),
!     &nresB(iB(m)),seqB(iB(m)),atomB(iB(m)),chB(iB(m)),
!     &xa(iA(m)),ya(iA(m)),za(iA(m)),xb(iB(m)),yb(iB(m)),zb(iB(m))
           endif
        enddo
      enddo
     
!      print*,'Number of Aligned Interfacial atoms:',in1
 
718   format(2(i4,1x,a3,1x,a3,1x,a1,5x),2x,3f8.3,2x,3f8.3)

!621   continue

      call u3b(w,r_1,r_2,in1,1,rms,u,t,ier) !u rotate r_1 to r_2
!      call u3b(w,r_1,r_2,n_ali,1,rms,u,t,ier) !u rotate r_1 to r_2
!      print*,'RMS: ',rms
      rmsd=dsqrt(rms/in1)
!      print*,'IRMS_progsub: ',rmsd,' for ',in1,' atoms'
      dt = u(1,1)*(u(2,2)*u(3,3) - u(2,3)*u(3,2))-u(1,2)*(u(2,1)*u(3,3)
     & - u(2,3)*u(3,1)) + u(1,3)*(u(2,1)*u(3,2) - u(2,2)*u(3,1))

      dmag = sqrt(t(1)**2 + t(2)**2 + t(3)**2)

!      write(*,88)'determinant of U : ',dt
!      print*,'magnitude of t : ',dmag

88    format(a25,f10.5)

!      write(*,211)u(1,1),u(1,2),u(1,3)
!      write(*,211)u(2,1),u(2,2),u(2,3)
!      write(*,211)u(3,1),u(3,2),u(3,3)
!      write(*,*)
!      write(*,211)t(1),t(2),t(3)
!      print*,'NUMBER OF INTERFACIAL ATOMS USED FOR SUPERPOSITION: ',in1
211   format(3(f12.3,2x))


!=============================================================================
            i1=2
!           write(*,1237)i1,atomB(i1),seqB(i1),chB(i1),nresB(i1),
!     &           xb(i1),yb(i1),zb(i1)

******************************************************************
******************************************************************
***   output rotated structure of pdb1:  superposing 1 onto 2 
***   output rotated structure of pdb1:  superposing 1 onto 2 
      sumsqd = 0.000
      Nt = 0
      idum = 0
!         OPEN(unit=7,file=outname1,status='unknown') !pdb1_tr
         do j=1,nseqA
            xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
            yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
            zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
!            write(7,1237)j,atomA(j),seqA(j),chA(j),nresA(j),
!     &           xt(j),yt(j),zt(j)
            iflag1 = 0
            iflag2 = 0
!================ GET THE MAPPED ATOM FROM THE NATIVE =============
            do i = 1,n_ali
                if (nresA(j)==nresB(iB(i)).and.chA(j).eq.chB(iB(i)))then
                if (atomA(j).eq.atomB(iB(i)))then
                i1 = iB(i)
                iflag1 = 1
                endif
                endif
            enddo
!=============== CHECK WHETHER THE ATOM BELONGS TO THE LIGAND ======
           if (iflag1 == 1)then
             do j1 = 1,iint
                if (nresB(i1)==nintAB(j1).and.chB(i1)==chAB(j1))then
                if (atomB(i1).eq.'N  '.or.atomB(i1).eq.'CA '.or.
     &atomB(i1).eq.'C  ')then
                iflag2 = 1
                endif
                endif
             enddo
          endif
          
!===================================================================
          if (iflag1 == 1 .and. iflag2 == 1)then
!           write(55,1237)i1,atomB(i1),seqB(i1),chB(i1),nresB(i1),
!     &           xb(i1),yb(i1),zb(i1)
!            write(88,1237)j,atomA(j),seqA(j),chA(j),nresA(j),
!     &           xt(j),yt(j),zt(j)
           d = sqrt((xt(j)-xb(i1))**2 + (yt(j)-yb(i1))**2
     &+ (zt(j)-zb(i1))**2)
!      write(224,561)nresB(i1),seqB(i1),chB(i1),atomB(i1),
!     &xb(i1),yb(i1),zb(i1),xt(j),yt(j),zt(j),d
             sumsqd = sumsqd + d**2
             Nt = Nt + 1
          endif
!====================================================================
         idum = idum + 1
         enddo
!         write(7,1238)
!         close(7)
 1237 format('ATOM  ',i5,2x,a3,1x,A3,1X,A1,1X,I3,4X,3F8.3)
 1238 format('TER')
!***   

!=======================================================================
!   RMSD OF THE ALIGNED INTERFACIAL ATOMS OF THE NATIVE AND THE 
!   TRANSFORMED MODEL
!=======================================================================

561   format(i4,2x,a3,1x,a1,2x,a3,2x,3f8.3,2x,3f8.3,2x,f10.3)

      irms = sqrt(sumsqd/float(Nt))

!      print*,'Nt= ',Nt,'IRMS= ',irms


!=============================================================================
!998   format('ONLY ',A1,' CHAIN CONSIDERED FOR SUPERPOSITION')

***   IRMSD calculation ---------------------------------->

!      print*,'LRMS:'
!      print*,chainsup,' IS TAKEN FOR SUPERPOSITION AND ',chaincal, 
!     &' FOR RMSD CALCULATION'

      do i = 1,3
      t(i) = 0.000
           do j = 1,3
           u(i,j) = 0.000
           enddo
      enddo

!================================================================

      in1 = 0

      do m=1,n_ali
        if ((chA(iA(m)).eq.chainsup).and.(chB(iB(m)).eq.chainsup))then
         in1 = in1 + 1
         r_1(1,in1)=xa(iA(m))
         r_1(2,in1)=ya(iA(m))
         r_1(3,in1)=za(iA(m))
         r_2(1,in1)=xb(iB(m))
         r_2(2,in1)=yb(iB(m))
         r_2(3,in1)=zb(iB(m))
!         print*,chA(iA(m)),'   ',chB(iB(m))
        endif
      enddo
      
!      print*,'Number of receptor atoms used for superposition: ',in1

!621   continue

!      call u3b(w,r_1,r_2,n_ali,1,rms,u,t,ier) !u rotate r_1 to r_2
      call u3b(w,r_1,r_2,in1,1,rms,u,t,ier) !u rotate r_1 to r_2
      rmsd=dsqrt(rms/in1)
!      print*,'Receptor-RMS_progsub: ',rmsd
      dt = u(1,1)*(u(2,2)*u(3,3) - u(2,3)*u(3,2))-u(1,2)*(u(2,1)*u(3,3)
     & - u(2,3)*u(3,1)) + u(1,3)*(u(2,1)*u(3,2) - u(2,2)*u(3,1))

      dmag = sqrt(t(1)**2 + t(2)**2 + t(3)**2)

******************************************************************
*     Output
******************************************************************
***   output RMSD ---------------------------->
******************************************************************
******************************************************************

******************************************************************
******************************************************************
      idum = 0
      sumsqd = 0.000
      Nt = 0

!         OPEN(unit=37,file=outname2,status='unknown') !pdb1_tr
         do j=1,nseqA
            xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
            yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
            zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
!            write(37,1237)j,atomA(j),seqA(j),chA(j),nresA(j),
!     &           xt(j),yt(j),zt(j)

            iflag1 = 0
            iflag2 = 0

!================ GET THE MAPPED ATOM FROM THE NATIVE =============
            do i = 1,n_ali
                if (nresA(j)==nresB(iB(i)).and.chA(j).eq.chB(iB(i)))then
                if (atomA(j).eq.atomB(iB(i)))then
                i1 = iB(i)
                iflag1 = 1
                endif
                endif
            enddo
!=============== CHECK WHETHER THE ATOM BELONGS TO THE LIGAND ======
                if (iflag1 == 1)then
                  if (chA(j).eq.chaincal)then
                  iflag2 = 1
                  endif
                endif
          if (iflag1 == 1 .and. iflag2 == 1)then
!           write(66,1237)i1,atomB(i1),seqB(i1),chB(i1),nresB(i1),
!     &           xb(i1),yb(i1),zb(i1)
!            write(99,1237)j,atomA(j),seqA(j),chA(j),nresA(j),
!     &           xt(j),yt(j),zt(j)
           d = sqrt((xt(j)-xb(i1))**2 + (yt(j)-yb(i1))**2
     &+ (zt(j)-zb(i1))**2)
!      write(324,561)nresB(i1),seqB(i1),chB(i1),atomB(i1),
!     &xb(i1),yb(i1),zb(i1),xt(j),yt(j),zt(j),d
             sumsqd = sumsqd + d**2
             Nt = Nt + 1
          endif
!====================================================================
         idum = idum + 1
         enddo
!         write(37,1238)
!         close(37)

126    format(3f8.3)
!       print*,idum

***   
!=====================================================================
!============= CALCULATE LRMS ========================================     
!============= RMS OVER THE PORTION THAT HAS NOT BEEN USED IN ========
!============= SUPERPOSITION =========================================
!=====================================================================

      lrms = sqrt(sumsqd/float(Nt))

!      print*,'Nt= ',Nt,'LRMS= ',lrms

!====================================================================
!     Calculate DockQ = avg(fnat, scaled_lrms, scaled_irms)
!====================================================================   

      d1 = 8.5
      d2 = 1.5

      scaled_lrms = 1/(1+(lrms/d1)**2)
      scaled_irms = 1/(1+(irms/d2)**2)
      DockQ = (fnat + scaled_lrms + scaled_irms)/float(3)

      c1 = 0.23
      c2 = 0.49
      c3 = 0.80

      capstat = ''
      print*,DockQ

      if(DockQ < c1)then
        capstat='Incorrect'
      elseif(DockQ >= c1 .and. DockQ < c2)then
        capstat='Acceptable'
          print*,'Here ?'
      elseif(DockQ >= c2 .and. DockQ < c3)then
        capstat='Medium'
      elseif(DockQ >= c3)then
        capstat='High'
      else
        capstat='Undef'
      endif

      print*,capstat

      write(245,234)fnat,lrms,irms,DockQ
      write(*,234)fnat,lrms,irms,DockQ
234   format('Fnat=',f8.3'  LRMS_bb=',f12.3,'  IRMS=',f12.3,
     &'  DockQ=',f12.5)
      write(245,246)capstat
      write(*,246)capstat

246   format('CAPRI_class (DockQ): ',a20)



!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================


!=====================================================================
!============= CALCULATE IRMS ========================================     
!============= RMS OVER THE INTERFACIAL BACKBONE ATOMS ALONE =========
!============= SUPERPOSE W.R.T. THE WHOLE MOLECULE ===================
!=====================================================================
!=====================================================================

******************************************************************
******************************************************************


***   output rotated structure of pdb1:  superposing 1 onto 2 
***   

555   continue

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END

      subroutine caldist(x1,y1,z1,x2,y2,z2,d)
      dxsq=(x1-x2)*(x1-x2)
      dysq=(y1-y2)*(y1-y2)
      dzsq=(z1-z2)*(z1-z2)
      d = sqrt(dxsq+dysq+dzsq)
      return
      end
     

ccccccccccc Calculate dRMSD between 1 and CA cccccccccccccccccccccccccccccccccc
      subroutine zydrmsd(r_1,r_2,nn,drms)
      double precision r_1(3,1),r_2(3,1),drms,d,dp
      drms=0
      do i=1,nn
         do j=i+1,nn
            d=sqrt((r_1(1,i)-r_1(1,j))**2+(r_1(2,i)-r_1(2,j))**2
     $           +(r_1(3,i)-r_1(3,j))**2)
            dp=sqrt((r_2(1,i)-r_2(1,j))**2+(r_2(2,i)-r_2(2,j))**2
     $           +(r_2(3,i)-r_2(3,j))**2)
            drms=drms+(d-dp)**2
         enddo
      enddo
      drms=sqrt(drms*2/(nn*(nn-1)))
      
      return
      end
      
cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
!      PARAMETER(nmax=100000)
!      double precision w(nmax), x(3,nmax), y(3,nmax)
      double precision w(*), x(3,*), y(3,*)
      integer n, mode
     
      double precision rms, u(3,3), t(3)
      integer ier
      
      integer i, j, k, l, m1, m
      integer ip(9), ip2312(4)
      double precision r(3,3), xc(3), yc(3), wc
      double precision a(3,3), b(3,3), e(3), rr(6), ss(6)
      double precision e0, d, spur, det, cof, h, g
      double precision cth, sth, sqrth, p, sigma
      
      double precision sqrt3, tol, zero
      
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
      
!      print*,'N from sub:',n

      wc = zero
      rms = zero
      e0 = zero
      
      do i=1, 3
         xc(i) = zero
         yc(i) = zero
         t(i) = zero
         do j=1, 3
            r(i,j) = zero
            u(i,j) = zero
            a(i,j) = zero
            if( i .eq. j ) then
               u(i,j) = 1.0
               a(i,j) = 1.0
            end if
         end do
      end do
      
      ier = -1
      if( n .lt. 1 ) return
      ier = -2
      do m=1, n
         if( w(m) .lt. 0.0 ) return
         wc = wc + w(m)
         do i=1, 3
            xc(i) = xc(i) + w(m)*x(i,m)
            yc(i) = yc(i) + w(m)*y(i,m)
         end do
      end do
      if( wc .le. zero ) return
      do i=1, 3
         xc(i) = xc(i) / wc
         yc(i) = yc(i) / wc
      end do
      
      do m=1, n
         do i=1, 3
            e0=e0+w(m)*((x(i,m)-xc(i))**2+(y(i,m)-yc(i))**2)
            d = w(m) * ( y(i,m) - yc(i) )
            do j=1, 3
               r(i,j) = r(i,j) + d*( x(j,m) - xc(j) )
            end do
         end do
      end do
      
      det = r(1,1) * ( (r(2,2)*r(3,3)) - (r(2,3)*r(3,2)) )
     &     - r(1,2) * ( (r(2,1)*r(3,3)) - (r(2,3)*r(3,1)) )
     &     + r(1,3) * ( (r(2,1)*r(3,2)) - (r(2,2)*r(3,1)) )
      
      sigma = det
      
      m = 0
      do j=1, 3
         do i=1, j
            m = m+1
            rr(m) = r(1,i)*r(1,j) + r(2,i)*r(2,j) + r(3,i)*r(3,j)
         end do
      end do
      
      spur = (rr(1)+rr(3)+rr(6)) / 3.0
      cof = (((((rr(3)*rr(6) - rr(5)*rr(5)) + rr(1)*rr(6))
     &     - rr(4)*rr(4)) + rr(1)*rr(3)) - rr(2)*rr(2)) / 3.0
      det = det*det
      
      do i=1, 3
         e(i) = spur
      end do
      if( spur .le. zero ) goto 40
      d = spur*spur
      h = d - cof
      g = (spur*cof - det)/2.0 - spur*h
      if( h .le. zero ) then
         if( mode .eq. 0 ) then
            goto 50
         else
            goto 30
         end if
      end if
      sqrth = dsqrt(h)
      d = h*h*h - g*g
      if( d .lt. zero ) d = zero
      d = datan2( dsqrt(d), -g ) / 3.0
      cth = sqrth * dcos(d)
      sth = sqrth*sqrt3*dsin(d)
      e(1) = (spur + cth) + cth
      e(2) = (spur - cth) + sth
      e(3) = (spur - cth) - sth
	
      if( mode .eq. 0 ) then
         goto 50
      end if
      
      do l=1, 3, 2
         d = e(l)
         ss(1) = (d-rr(3)) * (d-rr(6))  - rr(5)*rr(5)
         ss(2) = (d-rr(6)) * rr(2)      + rr(4)*rr(5)
         ss(3) = (d-rr(1)) * (d-rr(6))  - rr(4)*rr(4)
         ss(4) = (d-rr(3)) * rr(4)      + rr(2)*rr(5)
         ss(5) = (d-rr(1)) * rr(5)      + rr(2)*rr(4)
         ss(6) = (d-rr(1)) * (d-rr(3))  - rr(2)*rr(2)
         
         if( dabs(ss(1)) .ge. dabs(ss(3)) ) then
            j=1
            if( dabs(ss(1)) .lt. dabs(ss(6)) ) j = 3
         else if( dabs(ss(3)) .ge. dabs(ss(6)) ) then
            j = 2
         else
            j = 3
         end if
         
         d = zero
         j = 3 * (j - 1)
         
         do i=1, 3
            k = ip(i+j)
            a(i,l) = ss(k)
            d = d + ss(k)*ss(k)
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            a(i,l) = a(i,l) * d
         end do
      end do
      
      d = a(1,1)*a(1,3) + a(2,1)*a(2,3) + a(3,1)*a(3,3)
      if ((e(1) - e(2)) .gt. (e(2) - e(3))) then
         m1 = 3
         m = 1
      else
         m1 = 1
         m = 3
      endif
      
      p = zero
      do i=1, 3
         a(i,m1) = a(i,m1) - d*a(i,m)
         p = p + a(i,m1)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 21 i=1, 3
            if (p .lt. dabs(a(i,m))) goto 21
            p = dabs( a(i,m) )
            j = i
 21      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( a(k,m)**2 + a(l,m)**2 )
         if( p .le. tol ) goto 40
         a(j,m1) = zero
         a(k,m1) = -a(l,m)/p
         a(l,m1) =  a(k,m)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            a(i,m1) = a(i,m1)*p
         end do
      end if
      
      a(1,2) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
      a(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
      a(3,2) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
      
 30   do l=1, 2
         d = zero
         do i=1, 3
            b(i,l) = r(i,1)*a(1,l) + r(i,2)*a(2,l) + r(i,3)*a(3,l)
            d = d + b(i,l)**2
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            b(i,l) = b(i,l)*d
         end do
      end do
      d = b(1,1)*b(1,2) + b(2,1)*b(2,2) + b(3,1)*b(3,2)
      p = zero
      
      do i=1, 3
         b(i,2) = b(i,2) - d*b(i,1)
         p = p + b(i,2)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 22 i=1, 3
            if(p.lt.dabs(b(i,1)))goto 22
            p = dabs( b(i,1) )
            j = i
 22      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( b(k,1)**2 + b(l,1)**2 )
         if( p .le. tol ) goto 40
         b(j,2) = zero
         b(k,2) = -b(l,1)/p
         b(l,2) =  b(k,1)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            b(i,2) = b(i,2)*p
         end do
      end if
      
      b(1,3) = b(2,1)*b(3,2) - b(2,2)*b(3,1)
      b(2,3) = b(3,1)*b(1,2) - b(3,2)*b(1,1)
      b(3,3) = b(1,1)*b(2,2) - b(1,2)*b(2,1)
      
      do i=1, 3
         do j=1, 3
            u(i,j) = b(i,1)*a(j,1) + b(i,2)*a(j,2) + b(i,3)*a(j,3)
         end do
      end do
      
 40   do i=1, 3
         t(i) = ((yc(i) - u(i,1)*xc(1)) - u(i,2)*xc(2)) - u(i,3)*xc(3)
      end do
 50   do i=1, 3
         if( e(i) .lt. zero ) e(i) = zero
         e(i) = dsqrt( e(i) )
      end do
      
      ier = 0
      if( e(2) .le. (e(1) * 1.0d-05) ) ier = -1
      
      d = e(3)
      if( sigma .lt. 0.0 ) then
         d = - d
         if( (e(2) - e(3)) .le. (e(1) * 1.0d-05) ) ier = -1
      end if
      d = (d + e(2)) + e(1)
      
      rms = (e0 - d) - d
      if( rms .lt. 0.0 ) rms = 0.0
      
      return
      end



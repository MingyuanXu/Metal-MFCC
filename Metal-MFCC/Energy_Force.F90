subroutine Energy_Force_cal
use comparm
use mpi
implicit none

integer(kind=8)::i,j
real(kind=8)::tmp
if(.not.allocated(Force)) allocate(Force(3,natom))

open(1001,file='Force_information')
open(1002,file='Warning_information')
Energy=0
Force=0
call Double_counting
Energy=Energy+fenergy-cenergy-Dcount_E

write(101,'(A20,f20.8)')'1B_Energy:',Energy
do i=1,natom
    do j=1,3
        Force(j,i)=Force(j,i)+fforce(j,i)-cforce(j,i)-Dcount_F(j,i)
    end do
end do
write(1001,'(A)')'The Force array after the single body calculation'
do i=1,natom
   write(1001,'(i5,3f14.8)')i,(Force(j,i),j=1,3)
end do

if(center_atom.ne.0)then
    call QM_EF_correction
    Energy=Energy+QM_correct_E
    write(101,'(A20,f20.8)')'QM_correction:',QM_correct_E
    Force=Force+QM_correct_F
end if

if(center_atom.ne.0)then
    write(1001,'(A)')'the delta F will be added to Force array in subroutine of QM_correction'
    do i=1,natom
        write(1001,'(i5,3f14.8)')i,(QM_correct_F(j,i)*627.51/0.529,j=1,3)
    end do
    write(1001,'(A)')'The Force array after subroutine of QM_correction'
    do i=1,natom
        write(1001,'(i5,3f14.8)')i,(Force(j,i),j=1,3)
    end do
end if

call B2_EF_correction
Energy=Energy+Body2_Energy

write(101,'(A20,f20.8)')'2B_Energy:',Body2_Energy
if(b2_cutoff.eq.0)then
    call HBond_EF
end if
do i=1,natom
    do j=1,3
        Force(j,i)=Force(j,i)+Body2_Force(j,i)
    end do
end do
!D_pairs=D_pairs+B2D_pairs
do i=1,natom
    do j=i,natom
        if(D_pairs(i,j).ne.1)then
            write(1002,'(3i5)')i,j,D_pairs(i,j)
        end if
    end do
end do
write(1001,'(A)')'the delta F will be added to Force array in subroutine of B2_EF_correction'
do i=1,natom
    tmp=Body2_Force(1,i)**2+Body2_Force(2,i)**2+Body2_Force(3,i)**2
    tmp=sqrt(tmp)*627.51/0.529
    write(1001,'(i5,4f14.8)')i,(Body2_Force(j,i),j=1,3)
    if(tmp.gt.30)then
        write(1002,'(A16,i5,A22,f5.2,A14)')'delat F of atom ',i,' is too large,size is ',tmp,' kcal/mol/Ang!'
    end if
end do

write(1001,'(A)')'the Force array after subroutine of B2_EF_correction'
do i=1,natom
    write(1001,'(i5,3f14.8)')i,(Force(j,i),j=1,3)
end do
close(1001)
close(1002)
if(trim(adjustl(Frag_method)).ne.'EE-GMFCC')then
    call GMFCC_correct
else 
    do i=1,natom-1
        do j=i+1,natom
            if(GMFCC_pairs(i,j).eq.0)then
                tmpforce=0
                tmpenergy=0
                call cal_vdw(i,j)
                Energy=Energy+tmpenergy
                Force=Force+tmpforce
            end if
        end do
    end do
end if
if(add_bg.eq.'yes')then
    open(701,file='BG_Force.dat')
        do i=1,bg_num
            write(701,'(i10,3f14.8)')i,(bg_force(j,i),j=1,3)
        end do
    close(701)
end if
write(101,'(A16)')'Total energy is:'
write(101,'(f16.8)')Energy
tmp=0
write(*,'(A18,f14.8)')'MFCC join Forces: ',tmp
do i=1,natom
    do j=1,3
       tmp=tmp+Force(j,i)
    end do
end do
if(add_bg.eq.'yes')then
    do i=1,bg_num
        do j=1,3
            tmp=tmp+bg_force(j,i)
        end do
    end do
end if
tmp=tmp*627.51/0.529
write(*,'(A18,f20.8)')'MFCC join Forces: ',tmp
if(trim(adjustl(if_cal_charge)).eq.'yes')then
    write(101,'(A)')'FB MFCC charge'
    do i=1,natom
        write(101,'(i10,f14.8)')i,QM_charge(i)
    end do
end if
!Force=Force+restrain_F
end subroutine
!===================================================================================
subroutine AA_energy_force(i,j)
!this function is ok,nobug
use comparm
use mpi
implicit none

integer(kind=8)::i,j,k
real(kind=8)::dis
dis=atomdis(i,j)
dis=dis/rfactor
tmpenergy=0
tmpforce=0
tmpenergy=crg(i)*crg(j)/dis
do k=1,3
    tmpforce(k,i)=tmpforce(k,i)+crg(i)*crg(j)/dis**3*&
                    (atomcrd(k,i)-atomcrd(k,j))/rfactor
    tmpforce(k,j)=tmpforce(k,j)+crg(i)*crg(j)/dis**3*&
                    (atomcrd(k,j)-atomcrd(k,i))/rfactor
end do
if(coordinate_atom(i).eq.1)then
    do k=1,3
       tmpforce(k,i)=0
    end do
end if
if(coordinate_atom(j).eq.1)then
    do k=1,3
       tmpforce(k,j)=0
    end do
end if
end subroutine
!===================================================================================
subroutine judge_line_connect(i,j,line_flag)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,k,l
logical::line_flag
line_flag=.false.
do k=1,unit_num
    do l=1,unit_num
        if(QM_tag(k).eq.i.and.QM_tag(l).eq.j)then
            if(line_connect(k,l).eq.1.or.line_connect(l,k).eq.1)then
                line_flag=.true.
            end if
        end if
    end do
end do
end subroutine
!=================================================================================
subroutine GMFCC_correct
use comparm
use mpi
implicit none
integer(kind=8)::i,j
if(charge_type.ne.'Amoeba')then
    do i=1,natom
        do j=i+1,natom
            if(GMFCC_pairs(i,j).eq.0)then
                call AA_energy_force(i,j)
                call cal_vdw(i,j)
                Energy=Energy+tmpenergy
                Force=Force+tmpforce
            end if
        end do
    end do
else if(charge_type.eq.'Amoeba')then
    call print_GMFCC_pairs
    call cal_Amoeba_interaction
end if
end subroutine
!================================================================================

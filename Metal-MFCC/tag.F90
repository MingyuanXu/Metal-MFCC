subroutine make_QM_tag(i)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,k
do j=1,unit_num
    if(QM_tag(j).eq.i)then
        do k=unit_pt(1,j),unit_pt(2,j)
            gauss_tag(k)=1
        end do
    end if
end do
do j=1,unit_num-1
    if(QM_tag(j).eq.i.and.QM_tag(j+1).eq.i)then
        if(atomname(unit_pt(1,j)).eq.'C'.and.atomname(unit_pt(1,j+1)).eq.'CB')then
            do k=selectCA(resnum(unit_pt(1,j+1))),selectCB(resnum(unit_pt(1,j+1)))
                gauss_tag(k)=1
            end do
        else if(atomname(unit_pt(1,j)).eq.'CB'.and.atomname(unit_pt(1,j+1)).eq.'C')then
            do k=selectCA(resnum(unit_pt(1,j))),selectCB(resnum(unit_pt(1,j)))
                gauss_tag(k)=1
            end do
        end if
    end if
end do
end subroutine
!================================================================
subroutine make_gauss_tag(qmstart,qmend)

use comparm
use mpi
implicit none
integer(kind=8)::i,qmstart,qmend
gauss_tag=0
do i=qmstart,qmend
    if(unit_type(unit_index(i)).ne.'N')then
    gauss_tag(i)=1
    end if
end do

end subroutine
!===============================================================
subroutine make_u2_tag(i,j)

use comparm
use mpi
implicit none
integer(kind=8)::i,j,m

do m=unit_pt(1,i),unit_pt(2,i)
    gauss_tag(m)=1
end do
do m=unit_pt(1,j),unit_pt(2,j)
    gauss_tag(m)=1
end do
end subroutine
!==============================================================
subroutine make_b2_tag(i,j)

use comparm
use mpi
implicit none
integer(kind=8)::i,j,m
do m=1,natom
    if(fragment(unit_index(m)).eq.i.or.fragment(unit_index(m)).eq.j)then
        gauss_tag(m)=1
    end if
end do
end subroutine
!===============================================================
subroutine make_b1_tag(i)
use comparm
use mpi

implicit none
integer(kind=8)::i,m

do m=1,natom
    if(fragment(unit_index(m)).eq.i)then
        gauss_tag(m)=1
    end if
end do

end subroutine
!==============================================================
subroutine make_u1_tag(i)
use comparm
use mpi

implicit none
integer(kind=8)::i,m
do m=unit_pt(1,i),unit_pt(2,i)
    gauss_tag(m)=1
end do
end subroutine
!==============================================================
subroutine make_Recal_tag(i)
use comparm
use mpi
implicit none

integer(kind=8)::i,m!,n

do m=unit_pt(1,i),unit_pt(2,i)
    gauss_tag(m)=1
end do
!if(atomname(unit_pt(1,i)).eq.'CB')then
!    if(QM_tag(i+1).ne.0.or.QM_tag(i-1).ne.0)then
!        do n=selectCA(resnum(unit_pt(1,i))),unit_pt(1,i)
!           gauss_tag(n)=1
!        end do
!    end if
!end if
end subroutine

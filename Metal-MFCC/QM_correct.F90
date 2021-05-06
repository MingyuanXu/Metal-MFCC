subroutine QM_correct(job_type)
use comparm 
use mpi
implicit none
integer(kind=8)::i,j,m,n
character(len=1)::job_type
character(len=50)::filename
logical::tmp_tag(unit_num)

if(job_type.eq.'r')then
    if(.not.allocated(QM_1B_E))&
        allocate(QM_1B_E(QM_ligand))
        QM_1B_E=0
    if(.not.allocated(QM_1B_F))&
        allocate(QM_1B_F(3,natom,QM_ligand))
        QM_1B_F=0
    if(.not.allocated(QM_2B_E)) &
        allocate(QM_2B_E(QM_ligand*(QM_ligand-1)))
        QM_2B_E=0
    if(.not.allocated(QM_2B_F))&
        allocate(QM_2B_F(3,natom,QM_ligand*(QM_ligand-1)))
        QM_2B_F=0
    if(.not.allocated(QM_ALL_F)) &
        allocate(QM_ALL_F(3,natom))
        QM_ALL_F=0
        QM_ALL_E=0
end if

if(job_type.eq.'w')then
!unit_connect=line_connect
end if
do i=1,QM_ligand-1
    gauss_tag=0
    call make_QM_tag(i)
    write(filename,'(A10,I3.3)')'QM_ligand_',i
    if(job_type.eq.'w')then
        call g09_job(filename) 
    else if(job_type.eq.'r')then
        call read_data(filename)
        QM_1B_E(i)=tmpenergy-tmpcenergy
        QM_charge=QM_charge-tmpcharge
        do m=1,natom
            do n=1,3
                QM_1B_F(n,m,i)=tmpforce(n,m)
            end do
        end do
        if(add_bg.eq.'yes')then
            bg_force=bg_force-tmp_bgforce
        end if
    end if
end do
!if(strategy.eq.2)then
gauss_tag=0
do i=1,QM_ligand
    call make_QM_tag(i)
end do
write(filename,'(A9)')'QM_region'
if(job_type.eq.'w')then
    call g09_job(filename)
else if(job_type.eq.'r')then
    call read_data(filename)
    QM_ALL_E=tmpenergy-tmpcenergy
    QM_ALL_F=tmpforce
    QM_charge=QM_charge+tmpcharge
    if(add_bg.eq.'yes')then
        bg_force=bg_force+tmp_bgforce
    end if
end if
!end if
if(job_type.eq.'w')then
do i=1,unit_num-1
    do j=i+1,unit_num
        if(QM_tag(i).ne.0.and.QM_tag(j).ne.0.and.QM_tag(i).ne.QM_tag(j))then
            unit_connect(i,j)=unit_connect(i,j)+1
            unit_connect(j,i)=unit_connect(j,i)+1
        end if
    end do
end do
end if
if(job_type.eq.'w')then
    write(198,'(A)')'QM_correct'
    do i=1,unit_num
        write(198,'(80i2)')(unit_connect(j,i),j=1,unit_num)
    end do    
end if

!QM_2b_num=0
!do i=1,QM_ligand-2
!   do j=i+1,QM_ligand-1
!        QM_2b_num=QM_2b_num+1
!        gauss_tag=0
!        call make_QM_tag(i)
!        call make_QM_tag(j)
!        write(filename,'(A5,2i3.3)')'QM_2B',i,j
!        if(job_type.eq.'w')then
!           call g09_job(filename)
!        else if(job_type.eq.'r')then
!           call read_data(filename)
!           QM_2B_E(QM_2b_num)=tmpenergy-tmpcenergy
!           do m=1,natom
!               do n=1,3
!                  QM_2B_F(n,m,QM_2b_num)=QM_2B_F(n,m,QM_2b_num)+tmpforce(n,m)
!               end do
!           end do
!        end if
!        if(job_type.eq.'w')then
!            do m=1,unit_num
!                do n=1,unit_num
!                    if(QM_tag(m).eq.i.and.QM_tag(n).eq.j)then
!                        unit_connect(m,n)=unit_connect(m,n)-1
!                        unit_connect(n,m)=unit_connect(n,m)-1
!                    end if
!                end do
!            end do
!        end if
!   end do
!end do


!if(job_type.eq.'w')then
!    write(198,'(A)')'QM_cut2B'
!    do i=1,unit_num
!        write(198,'(80i2)')(unit_connect(j,i),j=1,unit_num)
!    end do    
!end if
tmp_tag=.false.
if(job_type.eq.'r')then
    if( .not.allocated(Re_2B_E)) allocate(Re_2B_E(Re_cal_num))
    if( .not.allocated(Re_2B_F)) allocate(Re_2B_F(3,natom,Re_cal_num))
    if( .not.allocated(Re_1B_E)) allocate(Re_1B_E(unit_num))
    if( .not.allocated(Re_1B_F)) allocate(Re_1B_F(3,natom,unit_num))
    Re_2B_E=0
    Re_2B_F=0
    Re_1B_E=0
    Re_1B_F=0
    if(.not.allocated(Re_connect)) allocate(Re_connect(2,Re_cal_num))
    Re_connect=0
end if
Re_cal_num=0
do i=1,unit_num-1
    do j=i+1,unit_num
        if(QM_tag(i).ne.0.and.QM_tag(j).ne.0.and.QM_tag(i).ne.QM_tag(j))then
            if(unit_connect(i,j).gt.1.or.unit_connect(j,i).gt.1)then
                gauss_tag=0
                Re_cal_num=Re_cal_num+1
                if(job_type.eq.'r')then
                    Re_connect(1,Re_cal_num)=i
                    Re_connect(2,Re_cal_num)=j
                end if
                tmp_tag(i)=.true.
                tmp_tag(j)=.true.
                call make_Recal_tag(i)
                call make_Recal_tag(j)
                write(filename,'(A6,i3.3,i3.3)')'Recal_',i,j
                if(job_type.eq.'w')then
                    call g09_job(filename)
                else if(job_type.eq.'r')then
                    call read_data(filename)
                    Re_2B_E(Re_cal_num)=tmpenergy-tmpcenergy
                    QM_charge=QM_charge-tmpcharge
                    do m=1,natom
                        do n=1,3
                            Re_2B_F(n,m,Re_cal_num)=tmpforce(n,m)
                        end do
                    end do
                    if(add_bg.eq.'yes')then
                    bg_force=bg_force-tmp_bgforce
                    end if
                end if
            end if
        end if
    end do
end do
do i=1,unit_num
    if(tmp_tag(i).eqv..true.)then
        gauss_tag=0
        call make_Recal_tag(i)
        write(filename,'(A6,i3.3)')'Recal_',i
        if(job_type.eq.'w')then
            call g09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            Re_1B_E(i)=tmpenergy-tmpcenergy
            QM_charge=QM_charge+tmpcharge
            do m=1,natom
                do n=1,3
                    Re_1B_F(n,m,i)=tmpforce(n,m)
                end do
            end do
            if(add_bg.eq.'yes')then
                bg_force=bg_force+tmp_bgforce
            end if
        end if
    end if
end do

end subroutine

subroutine QM_EF_correction
use comparm

implicit none

integer(kind=8)::i,j,k,m,n,tag1(natom),tag2(natom)
real(kind=8)::df(3),df1
if(.not.allocated(QM_correct_F)) allocate(QM_correct_F(3,natom))
QM_correct_E=0
QM_correct_F=0
do i=1,QM_ligand
    QM_ALL_E=QM_ALL_E-QM_1B_E(i)
    do m=1,natom
        do n=1,3
            QM_ALL_F(n,m)=QM_ALL_F(n,m)-QM_1B_F(n,m,i)
        end do
    end do
end do
QM_correct_E=QM_correct_E+QM_ALL_E
QM_correct_F=QM_correct_F+QM_ALL_F
do i=1,QM_ligand-1
    do j=i+1,QM_ligand
        gauss_tag=0
        call make_QM_tag(i)
        tag1=gauss_tag
        gauss_tag=0
        call make_QM_tag(j)
        tag2=gauss_tag
        do m=1,natom
            do n=1,natom
                if(tag1(m).ne.0.and.tag2(n).ne.0)then
                    if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                        call AA_energy_force(m,n)
                        QM_correct_E=QM_correct_E+tmpenergy
                        if(if_field.eq.'yes')then
                            QM_correct_F=QM_correct_F+tmpforce 
                        end if
                    end if
                end if
            end do
        end do
    end do
end do

!QM_2b_num=0
!do i=1,QM_ligand-2
!    do j=i+1,QM_ligand-1
!        QM_2b_num=QM_2b_num+1
!        QM_2B_E(QM_2b_num)=QM_2B_E(QM_2b_num)-QM_1B_E(i)-QM_1B_E(j)
!        do m=1,natom
!            do n=1,3
!                QM_2B_F(n,m,QM_2b_num)=QM_2B_F(n,m,QM_2b_num)-QM_1B_F(n,m,i)-QM_1B_F(n,m,j)
!            end do
!        end do
!        gauss_tag=0
!        call make_QM_tag(i)
!        tag1=gauss_tag
!        gauss_tag=0
!        call make_QM_tag(j)
!        tag2=gauss_tag
!        do m=1,natom
!            do n=1,natom
!                if(tag1(m).ne.0.and.tag2(n).ne.0)then
!                    call AA_energy_force(m,n)
!                    QM_2B_E(QM_2b_num)=QM_2B_E(QM_2b_num)+tmpenergy
!                    if(if_field.eq.'yes')then
!                        do k=1,natom
!                            do l=1,3
!                                QM_2B_F(l,k,QM_2b_num)=QM_2B_F(l,k,QM_2b_num)+tmpforce(l,k)
!                            end do
!                        end do         
!                    end if
!                end if
!            end do
!        end do
!!        QM_correct_E=QM_correct_E-QM_2B_E(QM_2b_num)
!        do m=1,natom
!            do n=1,3
!!                QM_correct_F(n,m)=QM_correct_F(n,m)-QM_2B_F(n,m,QM_2b_num)
!            end do
!        end do
!    end do
!end do

do i=1,Re_cal_num
    Re_2B_E(i)=Re_2B_E(i)-Re_1B_E(Re_connect(1,i))-Re_1B_E(Re_connect(2,i))
    QM_correct_E=QM_correct_E-Re_2B_E(i)
    do m=1,natom
        do n=1,3
            Re_2B_F(n,m,i)=Re_2B_F(n,m,i)-Re_1B_F(n,m,Re_connect(1,i))-&
                            Re_1B_F(n,m,Re_connect(2,i))
            QM_correct_F(n,m)=QM_correct_F(n,m)-Re_2B_F(n,m,i)
        end do
    end do
    df =0
    df1=0
    gauss_tag=0
    call make_Recal_tag(Re_connect(1,i))
    tag1=gauss_tag
    gauss_tag=0
    call make_Recal_tag(Re_connect(2,i))
    tag2=gauss_tag
    do m=1,natom
        if(tag1(m).ne.0)then
            do n=1,natom
                if(tag2(n).ne.0)then
                    if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                        call AA_energy_force(m,n)
                        QM_correct_E=QM_correct_E-tmpenergy
                        if(if_field.eq.'yes')then
                            QM_correct_F=QM_correct_F-tmpforce
                            do k=1,3
                                df(k)=df(k)+tmpforce(k,bad_point)
                            end do
                        end if
                    end if
                end if
            end do
        end if
    end do
    do k=1,3
        df(k)=df(k)+Re_2B_F(k,bad_point,i)
        df1=df1+df(k)**2
    end do
    df1=sqrt(df1)*627.51/0.529
    write(999,'(A8,3i3,f8.4)')'Recal_index',i,Re_connect(1,i),Re_connect(2,i),df1
end do
end subroutine

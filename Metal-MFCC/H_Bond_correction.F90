subroutine H_Bond_correction(job_type)
    use comparm
    use mpi
    implicit none
    integer(kind=8)::i,j,m,n
    logical::H_flag
    character(len=1)::job_type
    character(len=50)::filename

    if(job_type.eq.'w')then
        if(.not.allocated(HB_tag))      allocate(HB_tag(unit_num))
        HB_num=0
        HB_tag=.false.
        if(.not.allocated(HB_connect)) allocate(HB_connect(2,1000))
        do i=1,unit_num
            do j=i+1,unit_num
                if(unit_type(i).eq.'B'.and.unit_type(j).eq.'B')then
                     H_flag=.false.
                     call judge_HBond(i,j,H_flag)
                     if(H_flag.eqv..true.)then
                        HB_num=HB_num+1
                        HB_connect(1,HB_num)=i
                        HB_connect(2,HB_num)=j
                        HB_tag(i)=.true.
                        HB_tag(j)=.true.
                     end if
                end if
            end do
        end do
    end if

    if(job_type.eq.'r')then
        if(.not.allocated(HB_energy))     allocate(HB_energy(HB_num))
        if(.not.allocated(HB_Force))      allocate(HB_Force(3,natom,HB_num))
        if(.not.allocated(HB1_Force))     allocate(HB1_Force(3,natom,unit_num))
        if(.not.allocated(HB1_energy))    allocate(HB1_energy(unit_num))
        if(.not.allocated(HB_charge))    allocate(HB_charge(natom,HB_num))
        if(.not.allocated(HB1_charge))    allocate(HB1_charge(natom,unit_num))
        HB_energy=0
        HB_Force=0
        HB1_Force=0
        HB1_energy=0
        HB_charge=0
        HB1_charge=0
    end if    
    
    do i=1,HB_num
        gauss_tag=0
        call make_u2_tag(HB_connect(1,i),HB_connect(2,i))
        write(filename,'(A3,I3.3,A,I3.3,A2,i3.3)')'HB_',i,'R',HB_connect(1,i),'_R',HB_connect(2,i)
        if(job_type.eq.'w')then
            call g09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            HB_energy(i)=tmpenergy-tmpcenergy
            do m=1,natom
                do n=1,3
                    HB_Force(n,m,i)=HB_Force(n,m,i)+tmpforce(n,m)
                end do
                HB_charge(m,i)=HB_charge(m,i)+tmpcharge(m)
            end do
            if(add_bg.eq.'yes')then
                bg_force=bg_force+tmp_bgforce
            end if
        end if
    end do

    do i=1,unit_num
        if(HB_tag(i).eqv..true.)then
            gauss_tag=0
            call make_u1_tag(i)
            write(filename,'(A4,i3.3)')'HB1_',i
            if(job_type.eq.'w')then
                call g09_job(filename)
            else if(job_type.eq.'r')then
                call read_data(filename)
                HB1_energy(i)=tmpenergy-tmpcenergy
                do m=1,natom
                    do n=1,3
                        HB1_Force(n,m,i)=HB1_Force(n,m,i)+tmpforce(n,m)
                    end do
                HB1_charge(m,i)=HB1_charge(m,i)+tmpcharge(m)
                end do
                if(add_bg.eq.'yes')then
                    bg_force=bg_force-tmp_bgforce
                end if
            end if
        end if
    end do
end subroutine H_Bond_correction

subroutine HBond_EF
    use comparm
    use mpi
    implicit none
    integer(kind=8)::i,m,n
    real(kind=8)::E_HB
    real(kind=8),allocatable::F_HB(:,:)

    if(.not.allocated(F_HB))  allocate(F_HB(3,natom))

    do i=1,HB_num
        HB_energy(i)=HB_energy(i)-HB1_energy(HB_connect(1,i))&
            -HB1_energy(HB_connect(2,i))
        do m=1,natom
            do n=1,3
                HB_Force(n,m,i)=HB_Force(n,m,i)-HB1_Force(n,m,HB_connect(1,i))&
                    -HB1_Force(n,m,HB_connect(2,i))  
                F_HB(n,m)=F_HB(n,m)+HB_Force(n,m,i)
            end do
            HB_charge(m,i)=HB_charge(m,i)-HB1_charge(m,HB_connect(1,i))-HB1_charge(m,HB_connect(2,i))
            QM_charge(m)=QM_charge(m)+HB_charge(m,i)
        end do
        E_HB=E_HB+HB_energy(i)
        
    end do

    Energy=Energy+E_HB

    do m=1,natom
        do n=1,3
            Force(n,m)=Force(n,m)+F_HB(n,m)
        end do
    end do
end subroutine HBond_EF

subroutine judge_HBond(i,j,H_flag)
    use comparm
    use mpi
    implicit none
    integer(kind=8)::i,j,m,n
    real(kind=8)::dis,angle
    logical::H_flag
    H_flag=.false.
    if(line_connect(i,j).eq.0)then
        if(unit_pt(1,i).ne.center_atom.and.unit_pt(1,j).ne.center_atom)then
            if(frag_type(fragment(i)).eq.0.or.frag_type(fragment(j)).eq.0)then
                do m=unit_pt(1,i),unit_pt(2,i)
                    do n=unit_pt(1,j),unit_pt(2,j)
                        if(atomname(m).eq.'O'.and.atomname(n).eq.'H')then
                            dis=atomdis(m,n)
                            call cal_Angle(selectN(resnum(n)),n,m,angle)
                            if(dis.le.3.0.and.angle.gt.120)then
                                H_flag=.true.
                            end if
                        else if(atomname(m).eq.'H'.and.atomname(n).eq.'O')then
                            dis=atomdis(m,n)
                            call cal_Angle(selectN(resnum(m)),m,n,angle)
                            if(dis.le.3.0.and.angle.gt.120)then
                                H_flag=.true.
                            end if
                        end if
                    end do
                end do
            end if
        end if
    end if
end subroutine judge_HBond

subroutine Body2_interaction(job_type)
use comparm
use mpi

implicit none
integer(kind=8)::i,j,k,m,n
logical::b2_flag
character(len=1)::job_type
character(len=50)::filename
real(kind=8)::cutoff
if(job_type.eq.'w')then
    if(.not.allocated(B2_connect))   allocate(B2_connect(frag_num,frag_num))
    if(.not.allocated(B2_uconnect))  allocate(B2_uconnect(unit_num,unit_num))
    if(.not.allocated(B2_ufconnect)) allocate(B2_ufconnect(unit_num,frag_num))
    if(.not.allocated(B1_tag))       allocate(B1_tag(frag_num))
    if(.not.allocated(B1_utag))      allocate(B1_utag(unit_num))
    unit_connect=line_connect
    B2_uconnect=0
    B2_connect=0
    B2_ufconnect=0
    B1_utag=0
    B1_tag=0
    B2_unum=0
    B2_num=0
    B2_uf_num=0  
    do i=1,unit_num
        write(198,'(80i2)')(unit_connect(j,i),j=1,unit_num)
    end do
    
    do i=1,frag_num-1
        do j=i+1,frag_num
            b2_flag=.false.
            call judge_2b_connect(i,j,b2_flag)
            if(b2_flag.eqv..true.)then
                B2_num=B2_num+1
                B2_connect(i,j)=B2_num
                B2_connect(j,i)=B2_num
                B1_tag(i)=1
                B1_tag(j)=1
                do m=1,unit_num
                    do n=1,unit_num
                        if(fragment(m).eq.i.and.fragment(n).eq.j)then
                            unit_connect(m,n)=unit_connect(m,n)+1
                            unit_connect(n,m)=unit_connect(n,m)+1
                        end if
                    end do
                end do
            end if
        end do
    end do

    do i=1,unit_num
        write(198,'(80i2)')(unit_connect(j,i),j=1,unit_num)
    end do

    if(b2_cutoff.le.3.5)then
        cutoff=b2_cutoff
    else 
        cutoff=3.5
    end if

    do i=1,unit_num
        do j=1,frag_num
            b2_flag=.false.
            call judge_uf_connect(i,j,b2_cutoff,b2_flag)
            if(b2_flag.eqv..true.)then
                B2_uf_num=B2_uf_num+1
                B2_ufconnect(i,j)=B2_uf_num
                B1_tag(j)=1
                B1_utag(i)=1
                do k=1,unit_num
                    if(fragment(k).eq.j)then
                        unit_connect(i,k)=unit_connect(i,k)+1
                        unit_connect(k,i)=unit_connect(k,i)+1
                    end if
                end do
            end if
        end do
    end do

    do i=1,unit_num
        write(198,'(80i2)')(unit_connect(j,i),j=1,unit_num)
    end do

    if(b2_cutoff.le.3.5)then
        cutoff=b2_cutoff
    else 
        cutoff=3.5
    end if

    do i=1,unit_num-1
        do j=i+1,unit_num
            if(unit_connect(i,j).eq.0)then
                b2_flag=.false.
                call judge_u2_connect(i,j,cutoff,b2_flag)
                !if(i.eq.6 .and. j.eq.14)then
                !    b2_flag=.true.
                !end if
                if(b2_flag.eqv..true.)then
                    B2_unum=B2_unum+1
                    B2_uconnect(i,j)=B2_unum
                    B2_uconnect(j,i)=B2_unum
                    unit_connect(i,j)=unit_connect(i,j)+1
                    unit_connect(j,i)=unit_connect(j,i)+1
                    B1_utag(i)=1
                    B1_utag(j)=1
                end if
            end if
        end do
    end do

    do i=1,unit_num
        write(198,'(80i2)')(unit_connect(j,i),j=1,unit_num)
    end do

end if

if(job_type.eq.'r')then
    if(.not.allocated(B2_E))    allocate(B2_E(B2_num))
    if(.not.allocated(B2_F))    allocate(B2_F(3,natom,B2_num))
    if(.not.allocated(B1_E))    allocate(B1_E(unit_num))
    if(.not.allocated(B1_F))        allocate(B1_F(3,natom,frag_num))
    if(.not.allocated(B2_uE))       allocate(B2_uE(B2_unum))
    if(.not.allocated(B2_uF))       allocate(B2_uF(3,natom,B2_unum))
    if(.not.allocated(B1_uE))       allocate(B1_uE(unit_num))
    if(.not.allocated(B1_uF))       allocate(B1_uF(3,natom,unit_num))
    if(.not.allocated(B2_uf_E))     allocate(B2_uf_E(B2_uf_num))
    if(.not.allocated(B2_uf_F))     allocate(B2_uf_F(3,natom,B2_uf_num))
    if(.not.allocated(B2FD_pairs))  allocate(B2FD_pairs(natom,natom,B2_num))
    if(.not.allocated(B2D_pairs))   allocate(B2D_pairs(natom,natom))
    if(.not.allocated(B2UD_pairs))  allocate(B2UD_pairs(natom,natom,B2_unum))
    if(.not.allocated(B2UFD_pairs)) allocate(B2UFD_pairs(natom,natom,B2_uf_num))
    if(.not.allocated(B1FD_pairs))  allocate(B1FD_pairs(natom,natom,frag_num))
    if(.not.allocated(B1UD_pairs))  allocate(B1UD_pairs(natom,natom,unit_num))
    if(trim(adjustl(add_bg)).eq.'yes')then
        if(.not.allocated(B2F_BGF))      allocate(B2F_BGF(3,bg_num,B2_num))
        if(.not.allocated(B2U_BGF))      allocate(B2U_BGF(3,bg_num,B2_unum))
        if(.not.allocated(B2UF_BGF))      allocate(B2UF_BGF(3,bg_num,B2_uf_num))
        if(.not.allocated(B2_BGF))       allocate(B2_BGF(3,bg_num))
        if(.not.allocated(B1F_BGF))      allocate(B1F_BGF(3,bg_num,frag_num))
        if(.not.allocated(B1U_BGF))      allocate(B1U_BGF(3,bg_num,unit_num))
    end if
    B2_uE=0
    B2_uF=0
    B1_uE=0
    B1_uF=0
    B2_E=0
    B2_F=0
    B1_E=0
    B1_F=0
    B2_uf_E=0
    B2_uf_F=0

    B2FD_pairs=0
    B2D_pairs=0
    B2UD_pairs=0
    B2UFD_pairs=0
    B1FD_pairs=0
    B1UD_pairs=0
    B2F_BGF=0
    B2U_BGF=0
    B2UF_BGF=0
    B2_BGF=0
    B1F_BGF=0
    B1U_BGF=0

end if

do i=1,frag_num-1
    do j=i+1,frag_num
        if(B2_connect(i,j).ge.1)then
            gauss_tag=0
            call make_b2_tag(i,j)
            write(filename,'(A4,i3.3,A2,i3.3)')'B2_F',i,'_F',j
            if(job_type.eq.'w')then
                call g09_job(filename)
            else if(job_type.eq.'r')then
                call read_data(filename)
                QM_charge=QM_charge+tmpcharge
                B2_E(B2_connect(i,j))=tmpenergy-tmpcenergy
                do m=1,natom
                    do n=1,3
                        B2_F(n,m,B2_connect(i,j))=B2_F(n,m,B2_connect(i,j))+tmpforce(n,m)
                    end do
                end do
                if(trim(adjustl(add_bg)).eq.'yes')then
                do m=1,bg_num
                    do n=1,3
                        B2F_BGF(n,m,B2_connect(i,j))=B2F_BGF(n,m,B2_connect(i,j))+tmp_bgforce(n,m)
                    end do
                end do
                end if
                do m=1,natom
                    do n=1,natom
                        B2FD_pairs(m,n,B2_connect(i,j))=B2FD_pairs(m,n,B2_connect(i,j))+tmpD_pairs(m,n)
                    end do
                end do
            end if
        end if
    end do
end do

do i=1,unit_num
    do j=1,frag_num
        if(B2_ufconnect(i,j).ne.0)then
            gauss_tag=0
            call make_u1_tag(i)
            call make_b1_tag(j)
            write(filename,'(A4,i3.3,A2,i3.3)')'B2_U',i,'_F',j
            if(job_type.eq.'w')then
                call g09_job(filename)
            else if(job_type.eq.'r')then
                call read_data(filename)
                B2_uf_E(B2_ufconnect(i,j))=tmpenergy-tmpcenergy
                QM_charge=QM_charge+tmpcharge
                do m=1,natom
                    do n=1,3
                        B2_uf_F(n,m,B2_ufconnect(i,j))=tmpforce(n,m)+B2_uf_F(n,m,B2_ufconnect(i,j))
                    end do
                end do
                if(trim(adjustl(add_bg)).eq.'yes')then
                    do m=1,bg_num
                        do n=1,3
                            B2UF_BGF(n,m,B2_ufconnect(i,j))=B2UF_BGF(n,m,B2_ufconnect(i,j))+tmp_bgforce(n,m)
                        end do
                    end do
                end if
                do m=1,natom
                    do n=1,natom
                        B2UFD_pairs(m,n,B2_ufconnect(i,j))=B2UFD_pairs(m,n,B2_ufconnect(i,j))+tmpD_pairs(m,n)
                    end do
                end do
            end if
        end if
    end do
end do

do i=1,unit_num-1
    do j=i+1,unit_num
        if(B2_uconnect(i,j).ne.0)then
            gauss_tag=0
            call make_u2_tag(i,j)
            write(filename,'(A4,i3.3,a2,i3.3)')'B2_U',i,'_U',j
            if(job_type.eq.'w')then
                call g09_job(filename)
            else if(job_type.eq.'r')then
                call read_data(filename)
                B2_uE(B2_uconnect(i,j))=tmpenergy-tmpcenergy
                QM_charge=QM_charge+tmpcharge
                do m=1,natom
                    do n=1,3
                        B2_uF(n,m,B2_uconnect(i,j))=B2_uF(n,m,B2_uconnect(i,j))+tmpforce(n,m)
                    end do
                end do
                if(trim(adjustl(add_bg)).eq.'yes')then
                    do m=1,bg_num
                        do n=1,3
                            B2U_BGF(n,m,B2_uconnect(i,j))=B2U_BGF(n,m,B2_uconnect(i,j))+tmp_bgforce(n,m)
                        end do
                    end do
                end if
                do m=1,natom
                    do n=1,natom
                        B2UD_pairs(m,n,B2_uconnect(i,j))=B2UD_pairs(m,n,B2_uconnect(i,j))+tmpD_pairs(m,n)
                    end do
                end do
            end if
        end if
    end do
end do

do i=1,frag_num
    if(B1_tag(i).eq.1)then
        gauss_tag=0
        call make_b1_tag(i)
        write(filename,'(A4,i3.3)')'B1_F',i 
        if(job_type.eq.'w')then
            if(frag_type(i).ne.-1)then
            call g09_job(filename)
            end if
        else if(job_type.eq.'r')then
            if(frag_type(i).eq.-1)then
            write(filename,'(A8,i4.4)')'LIGFrag_',i
            end if
            call read_data(filename)
            B1_E(i)=tmpenergy-tmpcenergy
            QM_charge=QM_charge-tmpcharge
            do m=1,natom
                do n=1,3
                    B1_F(n,m,i)=B1_F(n,m,i)+tmpforce(n,m)
                end do
            end do
            do m=1,natom
                do n=1,natom
                    B1FD_pairs(m,n,i)=B1FD_pairs(m,n,i)+tmpD_pairs(m,n)
                end do
            end do
            if(trim(adjustl(add_bg)).eq.'yes')then
                do m=1,bg_num
                    do n=1,3
                        B1F_BGF(n,m,i)=B1F_BGF(n,m,i)+tmp_bgforce(n,m)
                    end do
                end do
            end if
        end if
    end if
end do

do i=1,unit_num
    if(B1_utag(i).eq.1)then
        gauss_tag=0
        call make_u1_tag(i)
        write(filename,'(A4,i3.3)')'B1_U',i
        if(job_type.eq.'w')then
            call g09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            B1_uE(i)=tmpenergy-tmpcenergy
            QM_charge=QM_charge-tmpcharge
            do m=1,natom
                do n=1,3
                    B1_uF(n,m,i)=B1_uF(n,m,i)+tmpforce(n,m)
                end do
            end do
            do m=1,natom
                do n=1,natom
                    B1UD_pairs(m,n,i)=B1UD_pairs(m,n,i)+tmpD_pairs(m,n)
                end do
            end do
            if(trim(adjustl(add_bg)).eq.'yes')then
                do m=1,bg_num
                    do n=1,3
                       B1U_BGF(n,m,i)=B1U_BGF(n,m,i)+tmp_bgforce(n,m)
                    end do
                end do
            end if
        end if
    end if
end do

end subroutine

subroutine B2_EF_correction
use comparm
use mpi
implicit none
integer(kind=8)::i,j,k,m,n
real(kind=8)::df(3),df1
write(999,*)'2B_correct'
if(.not.allocated(Body2_Force))  allocate(Body2_Force(3,natom))

Body2_Energy=0
Body2_Force=0
open(777,file='2B_correct_E')
write(777,'(A)')'F-F'
do i=1,frag_num-1
    do j=i+1,frag_num
        if(B2_connect(i,j).gt.0)then
            df=0
            df1=0
            write(777,'(2i10,3F14.8)')i,j,B2_E(B2_connect(i,j)),B1_E(i),B1_E(j)
            B2_E(B2_connect(i,j))=B2_E(B2_connect(i,j))-B1_E(i)-B1_E(j)
            write(777,'(2i10,3F14.8)')i,j,B2_E(B2_connect(i,j)),B1_E(i),B1_E(j)
            Body2_Energy=Body2_Energy+B2_E(B2_connect(i,j))
            do m=1,natom
                do n=1,3
                    B2_F(n,m,B2_connect(i,j))=B2_F(n,m,B2_connect(i,j))-B1_F(n,m,i)-B1_F(n,m,j)
                    Body2_Force(n,m)=Body2_Force(n,m)+B2_F(n,m,B2_connect(i,j))
                end do
                do n=1,natom
                    B2FD_pairs(m,n,B2_connect(i,j))=B2FD_pairs(m,n,B2_connect(i,j))&
                        -B1FD_pairs(m,n,i)-B1FD_pairs(m,n,j)
                    B2FD_pairs(n,m,B2_connect(i,j))=B2FD_pairs(n,m,B2_connect(i,j))&
                        -B1FD_pairs(n,m,i)-B1FD_pairs(n,m,j)
                    B2D_pairs(m,n)=B2D_pairs(m,n)+B2FD_pairs(m,n,B2_connect(i,j))
                    B2D_pairs(n,m)=B2D_pairs(n,m)+B2FD_pairs(n,m,B2_connect(i,j))
                end do
            end do
            if(trim(adjustl(add_bg)).eq.'yes')then
                do m=1,bg_num
                    do n=1,3
                        B2F_BGF(n,m,B2_connect(i,j))=B2F_BGF(n,m,B2_connect(i,j))-B1F_BGF(n,m,i)-B1F_BGF(n,m,j)
                        B2_BGF(n,m)=B2_BGF(n,m)+B2F_BGF(n,m,B2_connect(i,j))
                    end do
                end do
            end if
            do m=1,natom
                if(fragment(unit_index(m)).eq.i)then
                    do n=1,natom
                        if(fragment(unit_index(n)).eq.j)then
                            if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                            call AA_energy_force(m,n)
                               Body2_Energy=Body2_Energy+tmpenergy
                            if(if_field.eq.'yes')then
                                Body2_Force=Body2_Force+tmpforce
                                do k=1,3
                                    df(k)=df(k)+tmpforce(k,bad_point)
                                end do
                                B2D_pairs(m,n)=B2D_pairs(m,n)-1
                                B2D_pairs(n,m)=B2D_pairs(n,m)-1
                            end if
                            end if
                        end if
                    end do
                end if
            end do
            do n=1,3
                df(n)=B2_F(n,bad_point,B2_connect(i,j))+df(n)
                df1=df(n)**2+df1
            end do
            df1=sqrt(df1)*627.51/0.529
            write(999,'(A8,3i3,f8.4)')'B2_index',B2_connect(i,j),i,j,df1
        end if
    end do
end do
write(777,'(A)')'u-F'
do i=1,unit_num
    do j=1,frag_num
        if(B2_ufconnect(i,j).ne.0)then
            df=0
            df1=0
            write(777,'(2i10,3F14.8)')i,j,B2_uf_E(B2_ufconnect(i,j)),B1_uE(i),B1_E(j)
            B2_uf_E(B2_ufconnect(i,j))=B2_uf_E(B2_ufconnect(i,j))-B1_uE(i)-B1_E(j)
            Body2_Energy=Body2_Energy+B2_uf_E(B2_ufconnect(i,j))
            write(777,'(2i10,3F14.8)')i,j,B2_uf_E(B2_ufconnect(i,j)),B1_uE(i),B1_E(j)
            do m=1,natom
                do n=1,3
                    B2_uf_F(n,m,B2_ufconnect(i,j))=B2_uf_F(n,m,B2_ufconnect(i,j))-&
                                                   B1_uF(n,m,i)-B1_F(n,m,j)
                    Body2_Force(n,m)=Body2_Force(n,m)+B2_uf_F(n,m,B2_ufconnect(i,j))
                end do
                do n=1,natom
                    B2UFD_pairs(m,n,B2_ufconnect(i,j))=B2UFD_pairs(m,n,B2_ufconnect(i,j))&
                        -B1UD_pairs(m,n,i)-B1FD_pairs(m,n,j)
                    B2UFD_pairs(n,m,B2_ufconnect(i,j))=B2UFD_pairs(n,m,B2_ufconnect(i,j))&
                        -B1UD_pairs(n,m,i)-B1FD_pairs(n,m,j)
                    B2D_pairs(m,n)=B2D_pairs(m,n)+B2UFD_pairs(m,n,B2_ufconnect(i,j))
                    B2D_pairs(n,m)=B2D_pairs(n,m)+B2UFD_pairs(n,m,B2_ufconnect(i,j))
                end do
            end do
            if(trim(adjustl(add_bg)).eq.'yes')then
                do m=1,bg_num
                    do n=1,3
                        B2UF_BGF(n,m,B2_ufconnect(i,j))=B2UF_BGF(n,m,B2_ufconnect(i,j))-B1U_BGF(n,m,i)-B1F_BGF(n,m,j)
                        B2_BGF(n,m)=B2_BGF(n,m)+B2UF_BGF(n,m,B2_ufconnect(i,j))
                    end do
                end do
            end if
            do m=1,natom
                if(unit_index(m).eq.i)then
                    do n=1,natom
                        if(fragment(unit_index(n)).eq.j)then
                            if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                               call AA_energy_force(m,n)
                               Body2_Energy=Body2_Energy+tmpenergy
                               if(if_field.eq.'yes')then
                                   Body2_Force=Body2_Force+tmpforce
                                   do k=1,3
                                       df(k)=df(k)+tmpforce(k,bad_point)
                                   end do
                                   B2D_pairs(m,n)=B2D_pairs(m,n)-1
                                   B2D_pairs(n,m)=B2D_pairs(n,m)-1
                               end if
                            end if
                        end if
                    end do
                end if
            end do
            do n=1,3
                df(n)=B2_uf_F(n,bad_point,B2_ufconnect(i,j))+df(n)
                df1=df(n)**2+df1
            end do
            df1=sqrt(df1)*627.51/0.529
            write(999,'(A15,3i3,f8.4)')'B2_UF_index',B2_ufconnect(i,j),i,j,df1
        end if
    end do
end do
write(777,'(A)')'u-u'
do i=1,unit_num-1
    do j=i+1,unit_num
        if(B2_uconnect(i,j).gt.0)then
            df=0
            df1=0
            write(777,'(2i10,3F14.8)')i,j,B2_uE(B2_connect(i,j)),B1_uE(i),B1_uE(j)
            B2_uE(B2_uconnect(i,j))=B2_uE(B2_uconnect(i,j))-B1_uE(i)-B1_uE(j)
            Body2_Energy=Body2_Energy+B2_uE(B2_uconnect(i,j))
            write(777,'(2i10,3F14.8)')i,j,B2_uE(B2_connect(i,j)),B1_uE(i),B1_uE(j)
            do m=1,natom
                do n=1,3
                    B2_uF(n,m,B2_uconnect(i,j))=B2_uF(n,m,B2_uconnect(i,j))-B1_uF(n,m,i)-B1_uF(n,m,j)
                    Body2_Force(n,m)=Body2_Force(n,m)+B2_uF(n,m,B2_uconnect(i,j))
                end do
                do n=1,natom
                    B2UD_pairs(m,n,B2_uconnect(i,j))=B2UD_pairs(m,n,B2_uconnect(i,j))&
                        -B1UD_pairs(m,n,i)-B1UD_pairs(m,n,j)
                    B2UD_pairs(n,m,B2_uconnect(i,j))=B2UD_pairs(n,m,B2_uconnect(i,j))&
                        -B1UD_pairs(n,m,i)-B1UD_pairs(n,m,j)
                    B2D_pairs(m,n)=B2D_pairs(m,n)+B2UD_pairs(m,n,B2_uconnect(i,j))
                    B2D_pairs(n,m)=B2D_pairs(n,m)+B2UD_pairs(n,m,B2_uconnect(i,j))
                end do
            end do
            if(trim(adjustl(add_bg)).eq.'yes')then
                do m=1,bg_num
                    do n=1,3
                        B2U_BGF(n,m,B2_uconnect(i,j))=B2U_BGF(n,m,B2_uconnect(i,j))-B1U_BGF(n,m,i)-B1U_BGF(n,m,j)
                        B2_BGF(n,m)=B2_BGF(n,m)+B2U_BGF(n,m,B2_uconnect(i,j))
                    end do
                end do
            end if
            do m=unit_pt(1,i),unit_pt(2,i)
                do n=unit_pt(1,j),unit_pt(2,j)
                    if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                        call AA_energy_force(m,n)
                        Body2_Energy=Body2_Energy+tmpenergy
                        if(if_field.eq.'yes')then
                            Body2_Force=Body2_Force+tmpforce
                            do k=1,3
                                df(k)=df(k)+tmpforce(k,bad_point)
                            end do
                        end if
                        B2D_pairs(m,n)=B2D_pairs(m,n)-1
                        B2D_pairs(n,m)=B2D_pairs(n,m)-1
                    end if
                end do
            end do
            do n=1,3
                df(n)=B2_uF(n,bad_point,B2_uconnect(i,j))+df(n)
                df1=df(n)**2+df1
            end do
            df1=sqrt(df1)*627.51/0.529
            write(999,'(A8,3i3,f8.4)')'B2_uindex',B2_uconnect(i,j),i,j,df1
        end if
    end do
end do
end subroutine

subroutine judge_2B_connect(i,j,b2_flag)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,m,n
logical::b2_flag,flag
!real(kind=8)::cut_off
b2_flag=.false.
!if(frag_type(i).eq.1.or.frag_type(j).eq.1)then
!    cut_off=4.0
!else if(frag_type(i).eq.0.and.frag_type(j).eq.0)then
!    cut_off=b2_cutoff
!end if
if(frag_type(i).ne.1.or.frag_type(j).ne.1)then
    if(frag_type(i).ne.2.and.frag_type(j).ne.2)then
        do m=1,unit_num
            do n=1,unit_num
                if(fragment(m).eq.i.and.fragment(n).eq.j)then
                call judge_u2_connect(m,n,b2_cutoff,flag)
                if(flag.eqv..true.)then
                    b2_flag=.true.
                    write(789,'(4i5)')i,j,m,n
                end if
                end if
            end do
        end do
!        if(frag_type(i).eq.1.and.frag_type(j).eq.1)then
!            b2_flag=.true.
!        end if
        do m=1,unit_num
            do n=1,unit_num
                if(fragment(m).eq.i.and.fragment(n).eq.j)then
                    if(line_connect(m,n).eq.1)then
                        b2_flag=.false.
                    end if
                end if
            end do
        end do
        !For this spectial system -----start
!        if(i.eq.10 .and. j.eq.57)then
!            b2_flag=.false.
!        end if
        !------End
    end if
end if
end subroutine

subroutine judge_u2_connect(i,j,cutoff,b2_flag)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,m,n
logical::b2_flag
real(kind=8)::dis,cutoff

b2_flag=.false.

if(line_connect(i,j).eq.0)then
    if(unit_pt(1,i).ne.center_atom.and.unit_pt(1,j).ne.center_atom)then
        if(frag_type(fragment(i)).le.0.or.frag_type(fragment(j)).le.0)then
        do m=unit_pt(1,i),unit_pt(2,i)
            do n=unit_pt(1,j),unit_pt(2,j)
                dis=atomdis(m,n)
                if(dis.le.cutoff)then
                    b2_flag=.true.
                end if
            end do
        end do
        end if
    end if
end if

end subroutine

subroutine judge_uf_connect(i,j,cutoff,b2_flag)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,m
logical::b2_flag,flag
real(kind=8)::cutoff

b2_flag=.false.

if(frag_type(fragment(i)).le.0.or.frag_type(j).le.0)then
    do m=1,unit_num
        if(fragment(m).eq.j)then
            flag=.false.
            call judge_u2_connect(i,m,cutoff,flag)
            if(flag.eqv..true.)then
               b2_flag=.true.
            end if
        end if
    end do
    do m=1,unit_num
        if(fragment(m).eq.j)then
            if((line_connect(i,m).ne.0).or.(line_connect(m,i).ne.0))then
               b2_flag=.false.
            end if
            if(unit_connect(i,m).ne.0.or.unit_connect(m,i).ne.0)then
               b2_flag=.false.
            end if
        end if
    end do
end if

end subroutine

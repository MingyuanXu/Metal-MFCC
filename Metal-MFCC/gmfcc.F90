subroutine unitment_job(job_type)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,k,qmstart,qmend
character(len=50)::filename
character(len=1)::job_type
if(job_type.eq.'w')then
    if(.not.allocated(line_connect)) allocate(line_connect(unit_num,unit_num))
    line_connect=0
else if(job_type.eq.'r')then
    if(.not.allocated(fforce))  allocate(fforce(3,natom))
    if(.not.allocated(cforce))  allocate(cforce(3,natom))
    fforce=0
    cforce=0
end if
gauss_tag=0
do i=1,chain_num
    do j=chain_start(i)+1,chain_end(i)-1
        if(j.eq.(chain_start(i)+1))then
            qmstart = resstart(chain_start(i))
            qmend  = selectC(chain_start(i)+2)-1
        else if (j.eq.(chain_end(i)-1))then
            if(residue(j-1).eq.'PRO')then
                qmstart = selectN(j-1)
            else 
                qmstart = selectCA(j-1)
            end if
            qmend = resend(chain_end(i))
        else
            if(residue(j-1).eq.'PRO')then
                qmstart = selectN(j-1)
            else
                qmstart = selectCA(j-1)
            end if
            qmend = selectC(j+1)-1
        end if
        call make_gauss_tag(qmstart,qmend)
        call identify_unitment(qmstart,qmend)
        write(filename,'(A10,i3.3,a3,i3.3,a3)')&
                       'MFCC_chain',i,'_F_',j,residue(j)(1:3)
        if(job_type.eq.'w')then
            call G09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            fenergy=fenergy+tmpenergy-tmpcenergy
            fforce=fforce+tmpforce
            D_pairs=D_pairs+tmpD_pairs
            if(add_bg.eq.'yes')then
            bg_force=bg_force+tmp_bgforce
            end if
            QM_charge=QM_charge+tmpcharge
            write(123,'(A)')filename
            write(123,'(3f14.8)')tmpforce
        end if
    end do
    do j=chain_start(i)+1,chain_end(i)-2
        if(residue(j).eq.'PRO')then
            qmstart = selectN(j)
        else 
            qmstart = selectCA(j)
        end if
        qmend = selectC(j+1)-1
        call make_gauss_tag(qmstart,qmend)
        write(filename,'(A10,i3.3,a3,i3.3,a3)')&
                        'MFCC_chain',i,'_C_',j,residue(j)(1:3)
        if(job_type.eq.'w')then
            call G09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            cenergy=cenergy+tmpenergy-tmpcenergy
            cforce=cforce+tmpforce
            QM_charge=QM_charge-tmpcharge
            if(add_bg.eq.'yes')then
                bg_force=bg_force-tmp_bgforce
            end if
            D_pairs=D_pairs-tmpD_pairs
            write(123,'(A)')filename
            write(123,'(3f14.8)')tmpforce
        end if
    end do
end do
do i=1,combine_num
    gauss_tag=0
    do j=1,2
        qmstart=unit_pt(1,combine_unit(j,i))
        qmend=unit_pt(2,combine_unit(j,i))
        do k=qmstart,qmend
            gauss_tag(k)=1
        end do
    end do
    line_connect(combine_unit(1,i),combine_unit(2,i))=1
    line_connect(combine_unit(2,i),combine_unit(1,i))=1
    write(filename,'(A11,a3,I3.3,A1,a3,I3.3)')'SaltBridge_',&
                                    resname(unit_pt(1,combine_unit(1,i))),&
                                    resnum(unit_pt(1,combine_unit(1,i))),'_',&
                                    resname(unit_pt(1,combine_unit(2,i))),&
                                    resnum(unit_pt(1,combine_unit(2,i)))
    if(job_type.eq.'w')then
        call G09_job(filename)
    else if(job_type.eq.'r')then
        call read_data(filename)
        fenergy=fenergy+tmpenergy-tmpcenergy
        fforce=fforce+tmpforce
        QM_charge=QM_charge+tmpcharge
        if(add_bg.eq.'yes')then
            bg_force=bg_force+tmp_bgforce
        end if
        D_pairs=D_pairs+tmpD_pairs
        write(123,'(A)')filename
        write(123,'(3f14.8)')tmpforce
    end if
end do
do j=1,unit_num
    do k=1,unit_num
        if(fragment(j).eq.fragment(k).and.frag_type(fragment(j)).eq.-1)then
            line_connect(j,k)=1
            line_connect(k,j)=1
        end if
    end do
end do
do i=1,frag_num
    if(frag_type(i).eq.-1)then
        gauss_tag=0
        do j=1,natom
            if(fragment(unit_index(j)).eq.i)then
                gauss_tag(j)=1
            end if
        end do
        write(filename,'(A8,i4.4)')'LIGFrag_',i
        if(job_type.eq.'w')then
            call G09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            fenergy=fenergy+tmpenergy-tmpcenergy
            fforce=fforce+tmpforce
            QM_charge=QM_charge+tmpcharge
            if(add_bg.eq.'yes')then
                bg_force=bg_force+tmp_bgforce
            end if
        D_pairs=D_pairs+tmpD_pairs
        write(123,'(A)')filename
        write(123,'(3f14.8)')tmpforce
        end if
    end if
end do
!do i=1,nres
!    if(restype(i).eq.'L')then
!        qmstart=resstart(i)
!        qmend=resend(i)
!        call make_gauss_tag(qmstart,qmend)
!        write(filename,'(A11,i3.3,a3)')'MFCC_Ligand_',i,residue(i)(1:3)
!        if(job_type.eq.'w')then
!            call G09_job(filename)
!        else if(job_type.eq.'r')then
!            call read_data(filename)
!            fenergy=fenergy+tmpenergy-tmpcenergy
!            fforce=fforce+tmpforce
!        end if
!    end if
!end do
if(job_type.eq.'w')then
    unit_connect=line_connect
end if
return
end subroutine

subroutine unitment_job2(job_type)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,k,qmstart,qmend
character(len=50)::filename
character(len=1)::job_type
if(job_type.eq.'w')then
    if(.not.allocated(line_connect)) allocate(line_connect(unit_num,unit_num))
    line_connect=0
else if(job_type.eq.'r')then
    if(.not.allocated(fforce))  allocate(fforce(3,natom))
    if(.not.allocated(cforce))  allocate(cforce(3,natom))
    fforce=0
    cforce=0
end if
gauss_tag=0
do i=1,chain_num
    do j =chain_start(i),chain_end(i)
        if(j.eq.chain_start(i))then
            qmstart = resstart(chain_start(i))
            if(residue(j+1).eq.'PRO')then
                 qmend = selectC(j+1)-1
            else 
                 qmend = selectCA(j+1)-1
            end if
        else if(j.eq.chain_end(i))then
            qmend =resend(chain_end(i))
            qmstart =selectC(j-1)
        else 
            qmstart=selectC(j-1)
            if(residue(j+1).eq.'PRO')then
                qmend = selectC(j+1)-1
            else 
                qmend = selectN(j+1)+1
            end if
        end if
        call make_gauss_tag(qmstart,qmend)
        call identify_unitment(qmstart,qmend)
        write(filename,'(A10,i3.3,a3,i3.3,a3)')&
                       'MFCC_chain',i,'_F_',j,residue(j)(1:3)
        if(job_type.eq.'w')then
            call G09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            fenergy=fenergy+tmpenergy-tmpcenergy
            fforce=fforce+tmpforce
            QM_charge=QM_charge+tmpcharge
        end if
    end do
    do j=chain_start(i)+1,chain_end(i)
        qmstart = selectC(j-1)
        if(residue(j).eq.'PRO')then
            qmend=selectC(j)-1
        else 
            qmend=selectN(j)+1
        end if
        call make_gauss_tag(qmstart,qmend)
        write(filename,'(A10,i3.3,a3,i3.3,a3)')&
                    'MFCC_chain',i,'_C_',j,residue(j)(1:3)
        if(job_type.eq.'w')then
            call G09_job(filename)
        else if(job_type.eq.'r')then
            call read_data(filename)
            cenergy=cenergy+tmpenergy-tmpcenergy
            cforce=cforce+tmpforce
            QM_charge=QM_charge-tmpcharge
            write(123,'(A)')filename
            write(123,'(3f14.8)')tmpforce
        end if
    end do
end do
!do i=1,nres
!    if(restype(i).eq.'L')then
!        qmstart=resstart(i)
!        qmend=resend(i)
!        call make_gauss_tag(qmstart,qmend)
!        write(filename,'(A11,i3.3,a3)')'MFCC_Ligand_',i,residue(i)(1:3)
!        if(job_type.eq.'w')then
!            call G09_job(filename)
!        else if(job_type.eq.'r')then
!            call read_data(filename)
!            fenergy=fenergy+tmpenergy-tmpcenergy
!            fforce=fforce+tmpforce
!        end if
!    end if
!end do
do i=1,combine_num
    gauss_tag=0
    do j=1,2
        qmstart=unit_pt(1,combine_unit(1,i))
        qmend=unit_pt(2,combine_unit(1,i))
        do k=qmstart,qmend
            gauss_tag(k)=1
        end do
    end do
    line_connect(combine_unit(1,i),combine_unit(2,i))=1
    line_connect(combine_unit(2,i),combine_unit(1,i))=1
    
    write(filename,'(A11,a3,I3.3,A1,a3,I3.3)')'SaltBridge_',&
                                    resname(unit_pt(1,combine_unit(1,i))),&
                                    combine_unit(1,i),'_',&
                                    resname(unit_pt(1,combine_unit(2,i))),&
                                    combine_unit(2,i)
    if(job_type.eq.'w')then
        call G09_job(filename)
    else if(job_type.eq.'r')then
        call read_data(filename)
        fenergy=fenergy+tmpenergy-tmpcenergy
        fforce=fforce+tmpforce
        QM_charge=QM_charge+tmpcharge
        write(123,'(A)')filename
        write(123,'(3f14.8)')tmpforce
    end if
end do
if(job_type.eq.'w')then
    unit_connect=line_connect
end if
return
end subroutine

subroutine Double_counting
!this function is ok ,nobug
use comparm
use mpi
implicit none
integer(kind=8)::i,j,m,n,lmmstart,lmmend,rmmstart,rmmend

if(.not.allocated(Dcount_F))        allocate(Dcount_F(3,natom))
Dcount_E=0
Dcount_F=0
if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
!---------------------------------------------------
!one chain double counting
do i=1,chain_num
    do j=chain_start(i)+1,chain_end(i)-2
        if(j.eq.chain_start(i)+1)then
            lmmstart = resstart(chain_start(i))
        else 
            if(residue(j-1).eq.'PRO')then
                lmmstart = selectN(j-1)
            else
                lmmstart =selectCA(j-1)
            end if
        end if
        if(residue(j).eq.'PRO')then
            lmmend=selectC(j-1)+1
        else 
            lmmend=selectN(j)+1
        end if
        rmmstart=selectC(j+1)
        rmmend  =resend(chain_end(i))
        do m=lmmstart,lmmend
            do n=rmmstart,rmmend 
                if(unit_type(unit_index(m)).ne.'N'.and.unit_type(unit_index(m)).ne.'N')then
                    call AA_energy_force(m,n)
                    Dcount_E=Dcount_E+tmpenergy
                    if(trim(adjustl(if_field)).eq.'yes')then                
                            Dcount_F=Dcount_F+tmpforce
                    D_pairs(m,n)=D_pairs(m,n)-1
                    D_pairs(n,m)=D_pairs(n,m)-1
                    end if
                end if
            end do
        end do 
    end do
end do
!-------------------------------------------------
!chain-chain double counting
do i=1,chain_num-1
    do j=i+1,chain_num
        do m=resstart(chain_start(i)),resend(chain_end(i))
            do n=resstart(chain_start(j)),resend(chain_end(j))
                if(unit_type(unit_index(m)).ne.'N'.and.unit_type(unit_index(m)).ne.'N')then
                    call AA_energy_force(m,n)
                    Dcount_E=Dcount_E+tmpenergy
                    if(trim(adjustl(if_field)).eq.'yes')then                
                        Dcount_F=Dcount_F+tmpforce
                    D_pairs(m,n)=D_pairs(m,n)-1
                    D_pairs(n,m)=D_pairs(n,m)-1
                    end if
                end if
            end do
        end do
    end do
end do
!-------------------------------------------------
do i=1,natom-1
    if(unit_type(unit_index(i)).eq.'L'.or.unit_type(unit_index(i)).eq.'T')then
       do j=i+1,natom
            if(unit_type(unit_index(j)).ne.'L'.and.unit_type(unit_index(j)).ne.'T')then
                call AA_energy_force(i,j)
                Dcount_E=Dcount_E+tmpenergy
                if(trim(adjustl(if_field)).eq.'yes')then
                    Dcount_F=Dcount_F+tmpforce
                    D_pairs(i,j)=D_pairs(i,j)-1
                    D_pairs(j,i)=D_pairs(j,i)-1
                end if
            end if
            if(unit_type(unit_index(j)).eq.'L'.or.unit_type(unit_index(j)).eq.'T')then
                if(fragment(unit_index(i)).ne.fragment(unit_index(j)))then
                    call AA_energy_force(i,j)
                    Dcount_E=Dcount_E+tmpenergy
                    if(trim(adjustl(if_field)).eq.'yes')then
                        Dcount_F=Dcount_F+tmpforce
                    D_pairs(i,j)=D_pairs(i,j)-1
                    D_pairs(j,i)=D_pairs(j,i)-1
                    end if
                end if
            end if
       end do
    end if
end do
!=========================================================
do i=1,combine_num
    do m=1,natom
        if(unit_index(m).eq.combine_unit(1,i).or.&
            unit_index(m).eq.combine_unit(2,i))then
            if(unit_type(unit_index(m)).ne.'L'.and.unit_type(unit_index(m)).ne.'T')then
                do n=1,natom
                    if(unit_index(n).ne.combine_unit(1,i).and.&
                    unit_index(n).ne.combine_unit(2,i))then
                        if(unit_type(unit_index(n)).ne.'L'.and.unit_type(unit_index(n)).ne.'T')then
                            call AA_energy_force(m,n)
                            Dcount_E=Dcount_E+tmpenergy
                            if(trim(adjustl(if_field)).eq.'yes')then
                                Dcount_F=Dcount_F+tmpforce
                                D_pairs(m,n)=D_pairs(n,m)-1
                                D_pairs(n,m)=D_pairs(n,m)-1
                            end if
                        end if
                    end if
                end do
            end if
        end if
    end do
end do
end if
end subroutine

subroutine Orca_submit(filename)
      use comparm
      use mpi
      implicit none
      integer(kind=8)::i,j,charge,spin,H_num,electron_num,QM_atom
      character(len=30)::filename
      electron_num=0
      QM_atom=0
      do i=1,natom
          if(gauss_tag(i).eq.1)then
              QM_atom=QM_atom+1
          end if
      end do
      open(202,file=trim(adjustl(filename))//'.orcain')
      open(203,file=trim(adjustl(filename))//'.pdb')
          write(202,'(A10,i5)')'#QM_atom:',QM_atom
          write(202,'(A)')'!'//trim(adjustl(qmmethod))//' '//&
                trim(adjustl(qmbasis))//' '//trim(adjustl(add_word))
          if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
              open(204,file=trim(adjustl(filename))//'bg')
              write(202,'(A)')'%pointcharges'//&
                        trim(adjustl(filename))//'.bg'
          end if
          write(202,'(A)')'%pal nprocs '//trim(adjustl(qmproc))
          write(202,'(A)')'end'
          write(202,'(A)')'%maxcore '//trim(adjustl(qmmem))
          if(trim(adjustl(filename)).ne.'ALL_sys')then
              do i=1,natom
                 if(gauss_tag(i).eq.1)then
                     do j=i,natom
                         if(gauss_tag(j).eq.1)then
                             GMFCC_pairs(i,j)=1
                              GMFCC_pairs(j,i)=1
                          end if
                      end do
                  end if
              end do
          end if
          charge=0
          gauss_unit=0
          do i=1,natom
              if(gauss_tag(i).eq.1.and.unit_index(i).ne.0)then
                  gauss_unit(unit_index(i))=1
              end if
          end do
          do i=1,unit_num
              if(gauss_unit(i).eq.1)then
                  charge=charge+unit_charge(i)
              end if
          end do
          spin=1
          do i=1,natom
              if(gauss_tag(i).eq.1.and.elementnum(i).eq.29.and.&
                 resname(i).eq.'CU2')then
                 spin=2
             else if(gauss_tag(i).eq.1.and.elementnum(i).eq.26.and.&
                 resname(i).eq.'FE3')then
                 spin=2
             end if
          end do
          write(202,'(A6,i2,i2)')'* xyz ',charge,spin
          do i=1,natom
              if(gauss_tag(i).eq.1)then
                   write(202,'(A2,5x,3f14.8)')element(i),(atomcrd(j,i),j=1,3)
                   write(203,299) 'ATOM  ',i, atomname(i), resname(i)(1:3), &
                    resnum(i), (atomcrd(j,i),j=1,3)
                    electron_num=electron_num+elementnum(i)
              end if
          end do
          H_num=0
          do i=1,nbona
              if(gauss_tag(ib(i)).eq.1.and.gauss_tag(jb(i)).eq.0)then
                  call add_H(ib(i),jb(i))
                  H_num=H_num+1
                  electron_num=electron_num+1
              else if(gauss_tag(ib(i)).eq.0.and.gauss_tag(jb(i)).eq.1)then
                  call add_H(jb(i),ib(i))
                  H_num=H_num+1
                  electron_num=electron_num+1
              end if
          end do
          do i=1,nbonh
              if(gauss_tag(ibh(i)).eq.1.and.gauss_tag(jbh(i)).eq.0)then
                  call add_H(ibh(i),jbh(i))
                  H_num=H_num+1
                  electron_num=electron_num+1
              else if(gauss_tag(ibh(i)).eq.0.and.gauss_tag(jbh(i)).eq.1)then
                  call add_H(jbh(i),ibh(i))
                  H_num=H_num+1
                  electron_num=electron_num+1
              end if
          end do
          write(202,'(A)')'*'
          if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
              write(204,'(i5)')natom-QM_atom
              do i=1,natom
                  if(gauss_tag(i).eq.0)then
                      write(204,'(f14.8,5x,3f14.8)')crg(i),(atomcrd(j,i),j=1,3)
                  end if
              end do
              if(add_bg.eq.'yes')then
                  do i=1,bg_num
                      write(204,'(f14.8,5x,3f14.8)')bg_crg(i),(bg_crd(j,i),j=1,3)
                  end do
              end if
              close(204)
          end if
      close(202)
      close(203)
      if(trim(adjustl(filename)).ne.'ALL_sys')then
          write(102,'(A80,i20)')'orca  '//trim(adjustl(filename))//'.orcain >'&
              //trim(adjustl(filename))//'.out',electron_num**3
      end if
299 format(a6,1x,I4,1x,a4,1x,a3,2x,I4,4x,3f8.3)
end subroutine

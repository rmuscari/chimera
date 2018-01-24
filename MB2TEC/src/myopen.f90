!===============================================================================
subroutine myfopen(fileid,filename,fileform,filestatus)
implicit none

integer (kind=4), intent(in) :: fileid
character (len=*), intent(in) :: filename,fileform,filestatus

open (fileid,file=filename,form=fileform,status=filestatus,err=101)
rewind (fileid)

return

101 continue
write(*,'(/,a,a,a,/)') 'Errore nell''apertura del file "',trim(filename),'"'
stop

end

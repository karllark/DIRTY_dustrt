pro close_ps_output,ps=ps,eps=eps

if (keyword_set(ps) or keyword_set(eps)) then begin
    device,/close
    if (keyword_set(eps)) then begin
       if (eps NE 2) then begin
          set_plot,'x'
       endif
    endif else begin
       set_plot,'x'
    endelse
endif

end

;  14 Dec 2006 : added ability to ignore NaN and Infs
; ----------------------------------------------------------------------

function krange,karray,kplot_type=kplot_type

karray = double(karray)

if (not keyword_set(kplot_type)) then begin
   kplot_type = 'i'
endif

indxs = where(finite(karray),n_indxs)
if (n_indxs NE n_elements(karray)) then begin
    print,'Non-finite numbers in plotting array'
endif

if (n_indxs EQ 0) then begin
;    print,'No finite numbers in plotting array'
    return,[0.0,1.0]
endif else begin

    if (kplot_type EQ 'o') then begin
        indxs2 = where(karray(indxs) GT 0.0,n_indxs2)
        if (n_indxs2 GT 0) then begin
            indxs = indxs2[indxs]
        endif else begin
            print,'No array elements greater than 0'
        endelse
    endif
    
    kmm = dblarr(2)
    kmm(1) = max(karray[indxs],min=tflt)
    kmm(0) = tflt
    
    if (kplot_type EQ 'o') then begin
        kmm = alog10(kmm)
        width = kmm(1) - kmm(0)
        if (width EQ 0.0) then begin
            ret_vals = 10^[kmm[0] - 0.1,kmm[0] + 0.1]
        endif else begin
            ret_vals = [10^(kmm(0) - width*0.1),10^(kmm(1) + width*0.1)]
        endelse
        return,ret_vals
    endif else begin
        width = kmm(1) - kmm(0)
        if (width EQ 0.0) then begin
            if (kmm[0] NE 0.0) then begin
                ret_vals = kmm[0]*[0.9,1.1]
            endif else begin
                ret_vals = [-0.1,0.1]
            endelse
        endif else begin
            ret_vals = [kmm(0) - width*0.1,kmm(1) + width*0.1]
        endelse
        return,ret_vals
    endelse
endelse

end


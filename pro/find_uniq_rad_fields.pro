;+
; NAME:
;      FIND_UNIQ_RAD_FIELD
;
; PURPOSE:
;      Find the set of unique radiation fields
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      FIND_UNIQ_RAD_FIELD, filebase
;
; INPUTS:
;      filebase: filebase of FITS files output by DIRTYv2.
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; RESTRICTIONS:
;      This program written for research use.  No warrenties are
;      given.  Use at your own risk.
;
; PROCEDURE:
;     
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     Written by : Karl D. Gordon (4 Dec 2009)
;-
pro find_uniq_rad_fields, filebase,eps=eps,prompt=prompt

fits_read,filebase+'_rad_field.fits',rad_field,rad_header
size_cube = size(rad_field)
rad_field_unc = 0.25*rad_field

fits_read,filebase+'_wave_grid.fits',wave_grid,wave_header
wave_grid *= 1e4

n_min_points = fix(size_cube[4]/4.)

total_rad = fltarr(size_cube[1],size_cube[2],size_cube[3])

; loop over the radiation fields growing the set of unique fields
max_uniq = 500
uniq_indxs = intarr(3,max_uniq)
nonuniq_indxs = intarr(size_cube[1],size_cube[2],size_cube[3])
k_uniq = 0
same_chisqr = 3.
n_examined = 0
for i = 0,(size_cube[1]-1) do begin
    for j = 0,(size_cube[1]-1) do begin
        for k = 0,(size_cube[1]-1) do begin
            indxs = where(rad_field[i,j,k,*] GT 0,n_indxs)
            if (n_indxs GT n_min_points) then begin
                n_examined++
                chisqr = 1e6
                z = 0
                while ((chisqr GT same_chisqr) AND (z LT k_uniq)) do begin
                    indxs = where(rad_field[i,j,k,*]*rad_field[uniq_indxs[0,z],uniq_indxs[1,z],uniq_indxs[2,z],*] GT 0.0,n_indxs)
                    if (n_indxs GT 0) then begin
                        chisqr = total((rad_field[i,j,k,indxs] - rad_field[uniq_indxs[0,z],uniq_indxs[1,z],uniq_indxs[2,z],indxs])/ $
                                       rad_field_unc[uniq_indxs[0,z],uniq_indxs[1,z],uniq_indxs[2,z],indxs])^2
                        chisqr /= n_indxs
;                        print,z,i,j,k,chisqr
                    endif
                    z++
                endwhile
                if (chisqr GT same_chisqr) then begin ; new rad field
                    uniq_indxs[*,k_uniq] = [i,j,k]
                    nonuniq_indxs[i,j,k] = z
                    k_uniq++
                    if (k_uniq GE max_uniq) then begin
                        print,'more than '+strtrim(max_uniq,2)+' uniq radiation fields'
                        exit
                    endif
                endif else begin
                    nonuniq_indxs[i,j,k] = z - 1
                endelse
            endif else begin
                nonuniq_indxs[i,j,k] = -1
            endelse
        endfor
    endfor
endfor

print,'n unique rad fields = ' + strtrim(k_uniq,2) + '; n examined = ' + strtrim(n_examined,2)

; now plot the results

ans = ''
xrange = krange(wave_grid,kplot_type='o')
yrange = krange(rad_field,kplot_type='o')
yrange[1] = 1e2
yrange[0] = 1e-9

psym = 1
for i = 0,(k_uniq-1) do begin
    zindxs = where(nonuniq_indxs EQ i,n_zindxs)
    
    print,i,n_zindxs

    kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='oo', $
          xtitle='wavelength [micron]',ytitle='rad field'
    for k = 0,(n_zindxs-1) do begin
        mi = zindxs[k]
        save_mi = mi
        kk = mi / (size_cube[1]*size_cube[2])
        mi = mi mod (size_cube[1]*size_cube[2])
        jj = mi / size_cube[2]
        mi = mi mod size_cube[2]
        ii = mi
        indxs = where(rad_field[ii,jj,kk,*] GT 0,n_indxs)
        koplot,wave_grid[indxs],rad_field[ii,jj,kk,indxs],psym=psym,symsize=0.1

        indxs = where(rad_field[ii,jj,kk,*]*rad_field[uniq_indxs[0,i],uniq_indxs[1,i],uniq_indxs[2,i],*] GT 0.0,n_indxs)
        chisqr = total((rad_field[ii,jj,kk,indxs] - rad_field[uniq_indxs[0,i],uniq_indxs[1,i],uniq_indxs[2,i],indxs])/ $
                       rad_field_unc[uniq_indxs[0,i],uniq_indxs[1,i],uniq_indxs[2,i],indxs])^2
        chisqr /= n_indxs
        print,ii,jj,kk,uniq_indxs[0,i],uniq_indxs[1,i],uniq_indxs[2,i],k,chisqr
    endfor
    if (keyword_set(prompt)) then read,'Continue: ',ans
endfor

end

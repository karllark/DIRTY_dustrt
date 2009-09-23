;+
; NAME:
;      PLOT_RAD_FIELD_INFO
;
; PURPOSE:
;      Plot some information on the radiation field.
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      PLOT_RAD_FIELD_INFO, filebase
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
;     Written by : Karl D. Gordon (23 Sep 2009)
;-
pro plot_rad_field_info, filebase,eps=eps,prompt=prompt

fits_read,filebase+'_rad_field.fits',rad_field,rad_header
size_cube = size(rad_field)

fits_read,filebase+'_wave_grid.fits',wave_grid,wave_header
wave_grid *= 1e4

; now make the histogram
indxs = sort(abs(wave_grid - 0.55))
vindx = indxs[0]

; setup output
ps_file = filebase + '_norm_rad_field'
setup_ps_output,ps_file,eps=eps,ps=ps

; plot all the radiation fields on one plot
;window,xsize=800,ysize=600
xrange = krange(wave_grid,kplot_type='o')
yrange = krange(rad_field,kplot_type='o')
yrange[1] = 1e2
yrange[0] = 1e-9
kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='oo', $
      xtitle='wavelength [micron]',ytitle='rad field/mean(rad field)'
psym = 1
for i = 0,(size_cube[1]-1) do begin
    for j = 0,(size_cube[1]-1) do begin
        for k = 0,(size_cube[1]-1) do begin
            indxs = where(rad_field[i,j,k,*] GT 0,n_indxs)
            if ((n_indxs GT 0) AND (max(rad_field[i,j,k,*]) LT 1e10)) then begin
                mean_rad = mean(rad_field[i,j,k,indxs])
                koplot,wave_grid[indxs],rad_field[i,j,k,indxs]/mean_rad,psym=psym,symsize=0.1
            endif
        endfor
    endfor
endfor

close_ps_output,eps=eps,ps=ps

ans = ''
if (keyword_set(prompt)) then read,'Continue: ',ans

; setup output
ps_file = filebase + '_rad_field'
setup_ps_output,ps_file,eps=eps,ps=ps

; plot all the radiation fields on one plot
;window,xsize=800,ysize=600
xrange = krange(wave_grid,kplot_type='o')
yrange = krange(rad_field,kplot_type='o')
yrange[1] = 1e4
yrange[0] = 1e-9
kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='oo', $
      xtitle='wavelength [micron]',ytitle='rad field/mean(rad field)'
if (keyword_set(prompt)) then psym = 100 else psym = 1
for i = 0,(size_cube[1]-1) do begin
    for j = 0,(size_cube[1]-1) do begin
        for k = 0,(size_cube[1]-1) do begin
            indxs = where(rad_field[i,j,k,*] GT 0,n_indxs)
            if ((n_indxs GT 0) AND (max(rad_field[i,j,k,*]) LT 1e10)) then begin
                koplot,wave_grid[indxs],rad_field[i,j,k,indxs],psym=psym,symsize=0.1
                if (keyword_set(prompt)) then read,'Continue: ',ans
            endif
        endfor
    endfor
endfor

close_ps_output,eps=eps,ps=ps

if (keyword_set(prompt)) then read,'Continue: ',ans

vrad_field = rad_field[*,*,*,vindx]
good_indxs = where(vrad_field GT 0.0,n_good)
print,'n_good = ', n_good

min_val = min(vrad_field[good_indxs],max=max_val)
n_histo = 50
min_val_log = alog10(min_val)
max_val_log = alog10(max_val)
binsize = (max_val_log - min_val_log)/(n_histo)
histo = histogram(alog10(vrad_field[good_indxs]),min=min_val_log,max=max_val_log,binsize=binsize,locations=histo_vals,reverse_indices=r)
n_histo = n_elements(histo_vals)

histo_vals = 10^(histo_vals + (binsize/2.))

kplot,histo_vals,histo,kplot_type='oo'

; go through the bins and plot the rad fields
for i = 0,(n_histo-1) do begin
    print,i,histo[i]
    if (r[i] NE r[i+1]) then begin
        for k = 0,(histo[i]-1) do begin
            ; need to decode this: vrad_field[good_indxs[r[r[i]+k]]] 
            ; into i,j,k
        endfor
    endif

;    IF R[i] NE R[i+1] THEN A[R[R[I] : R[i+1]-1]] = 0 
endfor

end

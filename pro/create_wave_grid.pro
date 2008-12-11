; program to generate a wavelength grid
; written to allow a finer grid for the aromatic features 
; to be embedded in a coarser resolution full grid
;
; KDG - 9 Dec 2008

pro create_wave_grid,output_file=output_file

if (not keyword_set(output_file)) then output_file = 'wave_grid_pah_opt.dat'

coarse_min = 0.0912
coarse_max = 1e3
coarse_res = 10.

fine_min = 3.
fine_max = 20.
fine_res = 50.

n_coarse_waves = fix(alog10(coarse_max/coarse_min)/ $
                     alog10((1.0 + 2.0*coarse_res)/(2.0*coarse_res - 1.0)) + 1.0)
log_wave_min = alog10(coarse_min)
delta_log_wave = (alog10(coarse_max) - log_wave_min)/(n_coarse_waves - 1)

full_waves_log = log_wave_min + findgen(n_coarse_waves)*delta_log_wave
full_waves = 10^full_waves_log

; now get the fine grid
sindxs = sort(abs(full_waves - fine_min))
k1 = sindxs[0]
if (full_waves[k1] GT fine_min) then k1--
sindxs = sort(abs(full_waves - fine_max))
k2 = sindxs[0]
if (full_waves[k2] LT fine_max) then k2++

fine_min = full_waves[k1]
fine_max = full_waves[k2]

n_fine_waves = fix(alog10(fine_max/fine_min)/ $
                     alog10((1.0 + 2.0*fine_res)/(2.0*fine_res - 1.0)) + 1.0)
log_wave_min = alog10(fine_min)
delta_log_wave = (alog10(fine_max) - log_wave_min)/(n_fine_waves - 1)

fine_waves_log = log_wave_min + findgen(n_fine_waves)*delta_log_wave
fine_waves = 10^fine_waves_log

full_waves = [full_waves[0:k1-1],fine_waves,full_waves[k2+1:n_coarse_waves-1]]

openw,unit1,output_file,/get_lun
printf,unit1,'# wavelength grid optimized for pah features'
printf,unit1,'# created by create_wave_grid.pro'
printf,unit1,'#'
printf,unit1,'# micron'
n_waves = n_elements(full_waves)
for i = 0,(n_waves-1) do begin
    printf,unit1,full_waves[i]
endfor
free_lun,unit1

end

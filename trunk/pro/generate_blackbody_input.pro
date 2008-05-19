;+
; NAME:
;      GENERATE_BLACKBODY_INPUT
;
; PURPOSE:
;      Generate a stellar black body file for input to DIRTYv2.
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      GENERATE_BLACKBODY_INPUT
;
; INPUTS:
;      temperature : temperature in K
;      radius : radius in solar radii
;      out_filename : output filename
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;      ASCII file suitable for input to DIRTYv2
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
;     Written by : Karl D. Gordon (15 May 2008)
;-
pro generate_blackbody_input,temperature,radius,out_filename, $
  silent=silent

min_wave = 0.09
max_wave = 5000.0
wave_res = 100.

; get number of wavelength points
n_waves = fix(alog10(max_wave/min_wave)/ $
              alog10((1.0 + 2.0*wave_res)/(2.0*wave_res - 1.0)) + 1.0)
if (not keyword_set(silent)) then print,'n_waves = ',n_waves

; get the wavelength delta in log10 units
delta_log_wave = (alog10(max_wave) - alog10(min_wave))/(n_waves - 1)

; make the wavelength grid
waves = alog10(min_wave) + findgen(n_waves)*delta_log_wave
waves = 10^waves

if (not keyword_set(silent)) then begin
    print,'min wave = ',waves[0]
    print,'min wave = ',waves[n_waves-1]
endif

; get the blackbody
bb = blackbody_dirtyv2(waves,temperature)

; convert radius from solar radii
radius_m = radius*6.960d8

; convert to ergs s^-1 Hz^-1
bb_output = bb*1d7*4.d0*!PI*radius_m^2

; output file
openw,unit1,out_filename,/get_lun
printf,unit1,'#Blackbody from generate_blackbody_input.pro'
printf,unit1,'#temperature = ', temperature
printf,unit1,'#radius (solar radii) = ',  radius
printf,unit1,'#'
printf,unit1,'#wavelength [microns], luminosity [ergs s^-1 Hz^-1]'
for i = 0,(n_waves-1) do begin
    printf,unit1,waves[i],bb_output[i]
endfor
free_lun,unit1

end

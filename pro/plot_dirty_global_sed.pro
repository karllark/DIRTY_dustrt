;+
; NAME:
;      PLOT_DIRTY_GLOBAL_SED
;
; PURPOSE:
;      Plot the global SED(s) output by DIRTYv2 in an ASCII FITS table
;      file.
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      PLOT_DIRTY_GLOBAL_SED, filename
;
; INPUTS:
;      filename: name of ASCII FITS Table file output by DIRTYv2.
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
;     Written by : Karl D. Gordon (22 Mar 2008)
;-
pro plot_dirty_global_sed,filename,energy=energy,eps=eps, $
  xrange=xrange,yrange=yrange,table=table

table = mrdfits(filename,1,header)
table_tagnames = tag_names(table)

; print some statistics
tot_abs = fxpar(header,'TOT_ABS')
tot_emit = fxpar(header,'TOT_EMIT')
print,'total absorbed energy = ',tot_abs
print,'total emitted energy = ',tot_emit
print,'emit/abs = ',tot_emit/tot_abs

print,'# cells with enough = ',fxpar(header,'N_ENOUGH')
print,'# cells with not enough = ',fxpar(header,'N_NOT_EN')
print,'# cells with zero = ',fxpar(header,'N_ZERO')
print,'# cells with too few waves = ',fxpar(header,'N_FEWWAV')

xtitle = '!4k!3 [!4l!3m]'
if (keyword_set(energy)) then begin
    ytitle = '!4k!3F(!4k!3) [ergs s!U-1!N]'
    table.flux *= table.wavelength
    table.flux_input *= table.wavelength
    table.flux_rt_d *= table.wavelength
    table.flux_rt_s *= table.wavelength
    sindxs = where(table_tagnames EQ 'FLUX_ERE_D',n_sindxs)
    if (n_sindxs GT 0) then begin
        table.flux_ere_d *= table.wavelength
        table.flux_ere_s *= table.wavelength
    endif
    sindxs = where(table_tagnames EQ 'FLUX_DE_D_1',n_sindxs)
    if (n_sindxs GT 0) then begin
        table.flux_de_d_1 *= table.wavelength
        table.flux_de_s_1 *= table.wavelength
        table.flux_de_d_2 *= table.wavelength
        table.flux_de_s_2 *= table.wavelength
        table.flux_de_d_3 *= table.wavelength
        table.flux_de_s_3 *= table.wavelength
        
        sindxs = where(table_tagnames EQ 'FLUX_DE_D_4',n_sindxs)
        if (n_sindxs GT 0) then begin
            table.flux_de_d_4 *= table.wavelength
            table.flux_de_s_4 *= table.wavelength
            table.flux_de_d_5 *= table.wavelength
            table.flux_de_s_5 *= table.wavelength
        endif
        sindxs = where(table_tagnames EQ 'FLUX_DE_D_6',n_sindxs)
        if (n_sindxs GT 0) then begin
            table.flux_de_d_6 *= table.wavelength
            table.flux_de_s_6 *= table.wavelength
            table.flux_de_d_7 *= table.wavelength
            table.flux_de_s_7 *= table.wavelength
        endif
    endif
endif else begin
    ytitle = 'F(!4k!3) [ergs s!U-1!N !4l!3m!U-1!N]'
endelse

if (not keyword_set(xrange)) then xrange = krange(table.wavelength,kplot_type='o')
if (not keyword_set(yrange)) then yrange = krange([table.flux,table.flux_input],kplot_type='o')

setup_ps_output,repstr(filename,'.fits',''),eps=eps,bw=bw

setup_colors,base_color,back_color,blue_color,red_color,green_color, $
             yel_color,purple_color,light_blue_color,grey_color, $
             line_color=line_color,red=red,green=green,blue=blue, $
             bw=bw

; setup the plot
kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='oo', $
      xtitle=xtitle,ytitle=ytitle,color=base_color,background=back_color

; plot the input SED
koplot,table.wavelength,table.flux_input,psym=100,color=base_color,linestyle=2

; plot the RT components
koplot,table.wavelength,table.flux_rt_d+table.flux_rt_s,psym=100,color=blue_color,linestyle=0
koplot,table.wavelength,table.flux_rt_d,psym=100,color=blue_color,linestyle=1
koplot,table.wavelength,table.flux_rt_s,psym=100,color=blue_color,linestyle=2

sindxs = where(table_tagnames EQ 'FLUX_ERE_D',n_sindxs)
if (n_sindxs GT 0) then begin
    koplot,table.wavelength,table.flux_ere_d+table.flux_ere_s,psym=100,color=red_color,linestyle=0
    koplot,table.wavelength,table.flux_ere_d,psym=100,color=red_color,linestyle=1
    koplot,table.wavelength,table.flux_ere_s,psym=100,color=red_color,linestyle=2
endif

sindxs = where(table_tagnames EQ 'FLUX_DE_D_1',n_sindxs)
if (n_sindxs GT 0) then begin
; plot the total DE
    koplot,table.wavelength,table.flux_de_d_1+table.flux_de_s_1,psym=100,color=red_color,linestyle=0
    koplot,table.wavelength,table.flux_de_d_1,psym=100,color=red_color,linestyle=1
    koplot,table.wavelength,table.flux_de_s_1,psym=100,color=red_color,linestyle=2

; plot the DE component 
    koplot,table.wavelength,table.flux_de_d_2,psym=100,color=green_color,linestyle=1
    koplot,table.wavelength,table.flux_de_s_2,psym=100,color=green_color,linestyle=2
    koplot,table.wavelength,table.flux_de_d_3,psym=100,color=green_color,linestyle=3
    koplot,table.wavelength,table.flux_de_s_3,psym=100,color=green_color,linestyle=4
endif

sindxs = where(table_tagnames EQ 'FLUX_DE_D_4',n_sindxs)
if (n_sindxs GT 0) then begin
; plot the DE component 
    koplot,table.wavelength,table.flux_de_d_4,psym=100,color=light_blue_color,linestyle=1
    koplot,table.wavelength,table.flux_de_s_4,psym=100,color=light_blue_color,linestyle=2
    koplot,table.wavelength,table.flux_de_d_5,psym=100,color=light_blue_color,linestyle=3
    koplot,table.wavelength,table.flux_de_s_5,psym=100,color=light_blue_color,linestyle=4
endif

sindxs = where(table_tagnames EQ 'FLUX_DE_D_6',n_sindxs)
if (n_sindxs GT 0) then begin
; plot the DE component 
;purple_color=base_color
    koplot,table.wavelength,table.flux_de_d_6,psym=100,color=purple_color,linestyle=1
    koplot,table.wavelength,table.flux_de_s_6,psym=100,color=purple_color,linestyle=2
    koplot,table.wavelength,table.flux_de_d_7,psym=100,color=purple_color,linestyle=3
    koplot,table.wavelength,table.flux_de_s_7,psym=100,color=purple_color,linestyle=4
endif

; plot the total SED
koplot,table.wavelength,table.flux,psym=100,color=base_color

close_ps_output,eps=eps

end

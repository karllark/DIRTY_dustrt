;+
; NAME:
;      COMP_DIRTY_GLOBAL_SEDS
;
; PURPOSE:
;      Copmare the global SED(s) output by multiple runs of the 
;      DIRTYv2 model.
;
; CATEGORY:
;      DIRTYv2.
;
; CALLING SEQUENCE:
;      COMP_DIRTY_GLOBAL_SEDS, filenames
;
; INPUTS:
;      filenames: vector of the names of ASCII FITS Table file output by DIRTYv2.
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
;     Written by : Karl D. Gordon (25 Aug 2008)
;-
pro comp_dirty_global_seds,filenames,energy=energy,eps=eps, $
  xrange=xrange,yrange=yrange

n_files = n_elements(filenames)

setup_ps_output,repstr(filenames[0],'.fits','_comps'),eps=eps,bw=bw
    
for i = 0,(n_files-1) do begin
    table = mrdfits(filenames[i],1,header)
    table_tagnames = tag_names(table)

    xtitle = '!4k!3 [!4l!3m]'
    if (keyword_set(energy)) then begin
        ytitle = '!4k!3F(!4k!3) [!4l!3m!U-1!N ergs s!U-1!N !4l!3m!U-1!N]'
        table.flux *= table.wavelength
        table.flux_input *= table.wavelength
        table.flux_rt_d *= table.wavelength
        table.flux_rt_s *= table.wavelength
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

    if (i EQ 0) then begin
    ; setup the plot
        setup_colors,base_color,back_color,blue_color,red_color,green_color, $
                     yel_color,purple_color,light_blue_color,grey_color, $
                     line_color=line_color,red=red,green=green,blue=blue, $
                     bw=bw

        kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,kplot_type='oo', $
              xtitle=xtitle,ytitle=ytitle,color=base_color,background=back_color
    endif

    ; plot the total SED
    koplot,table.wavelength,table.flux,psym=100,color=line_color[i mod 5]

endfor

close_ps_output,eps=eps

end

pro setup_ps_output,output_filename,ps=ps,eps=eps,bw=bw, $
                    square=square

if (keyword_set(ps) or keyword_set(eps)) then begin
    set_plot,'ps'
    if (keyword_set(bw)) then begin
        output_filename = output_filename + '_bw'
    endif
    if (keyword_set(ps)) then begin
        device,filename=output_filename+'.ps'
        device,/landscape
        device,encapsulated=0
    endif else begin
        device,filename=output_filename+'.eps'
        device,/portrait
        device,encapsulated=1
    endelse
    device,/color
    if (keyword_set(square)) then begin
        device,xsize=10,ysize=10,/inches
    endif
    device,bits_per_pixel=8
endif

end

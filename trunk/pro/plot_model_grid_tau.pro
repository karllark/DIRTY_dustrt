
pro plot_model_grid_tau,filebase,kxrange=kxrange,kyrange=kyrange

fits_open,filebase+'_pos.fits',fcb_pos
fits_open,filebase+'_tau_ref_per_pc.fits',fcb_tau

for m = 1,fcb_pos.nextend do begin
    fits_read,fcb_pos,pos,header_pos,exten_no=m
    fits_read,fcb_tau,tau,header_tau,exten_no=m

    size_cube = size(tau)
    ; generate the radius cube
    radius = fltarr(size_cube[1],size_cube[2],size_cube[3])
    for i = 0,(size_cube[1]-1) do begin
        x_val = (pos[0,i] + pos[0,i+1])/2.0;
        for j = 0,(size_cube[2]-1) do begin
            y_val = (pos[1,j] + pos[1,j+1])/2.0 ;
            for k = 0,(size_cube[3]-1) do begin
                z_val = (pos[2,k] + pos[2,k+1])/2.0 ;
                radius[i,j,k] = sqrt(x_val^2 + y_val^2 + z_val^2)
            endfor
        endfor
    endfor

    if (m EQ 1) then begin
        if (not keyword_set(kxrange)) then kxrange = krange(radius,kplot_type='o')
        if (not keyword_set(kyrange)) then kyrange = krange(tau,kplot_type='o')
        kplot,[1],[1],/no_data,xrange=kxrange,yrange=kyrange,kplot_type='oo'
    endif
    
    koplot,radius,tau,psym=1
endfor

fits_close,fcb_pos
fits_close,fcb_tau

end

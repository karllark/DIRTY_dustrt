
pro setup_colors,base_color,back_color,blue_color,red_color,green_color, $
                 yel_color,purple_color,light_blue_color,grey_color, $
                 line_color=line_color,red=red,green=green,blue=blue, $
                 bw=bw

if (not keyword_set(bw)) then begin
    base_color = 1
    back_color = 0
    blue_color = 4
    red_color = 2
    green_color = 3
    yel_color = 5
    purple_color = 6
    light_blue_color = 7
    grey_color = 9
    line_color = [green_color,blue_color,red_color,purple_color, $
                  light_blue_color]
    
    red =   [1,0,1,0,0,0,1,0,0,0.75]
    green = [1,0,0,1,0,0.5,0,1,0.5,0.75]
    blue =  [1,0,0,0,1,0,1,1,0,0.75]
    tvlct,255*red,255*green,255*blue
endif else begin
    loadct,0,/silent
    base_color = 1
    back_color = 255
    blue_color = base_color
    red_color = base_color
    green_color = base_color
    yel_color = base_color
    purple_color = base_color
    light_blue_color = base_color
    line_color = replicate(base_color,10)
endelse

end

; ----------------------------------------------------------------------

pro get_new_krange,in_krange,in_karray,kplot_type=kplot_type

if (not keyword_set(kplot_type)) then begin
   kplot_type = 'i'
endif

krange2 = krange(in_karray,kplot_type=kplot_type)

new_krange = in_krange
new_krange(0) = min([in_krange(0),krange2(0)])
new_krange(1) = max([in_krange(1),krange2(1)])
in_krange = new_krange

end


a=''
start:
l=findfile("swp*mxhi")
readmx,l(0),h,wave,flux,wrange=[1350,1450]
jd_start=sxpar(h,'SJD-OBS')
plot,wave,smooth(flux,3),yr=[0,4e-8]
for i=0,n_elements(l)-1 do begin          
 readmx,l(i),h2,wave2,flux2,wrange=[1350,1450]
 jd_start2=sxpar(h2,'SJD-OBS')
 plot,wave,smooth(flux,3),yr=[0,4e-8],title=strtrim(jd_start,2)+" vs "+strtrim(jd_start2,2)
 oplot,wave2,smooth(flux2,3),col=220
 wait,1                            
 endfor 
goto,start
end

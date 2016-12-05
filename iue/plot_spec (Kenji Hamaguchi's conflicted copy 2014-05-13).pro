a=''
start:
l=findfile("swp*mxhi")
readmx,l(0),h,wave,flux,wrange=[1350,1450]
plot,wave,smooth(flux,3),yr=[0,4e-8]
for i=0,n_elements(l)-1 do begin          
 plot,wave,smooth(flux,3),yr=[0,4e-8],title=l(i)
 readmx,l(i),h,wave,flux,wrange=[1350,1450]
 oplot,wave,smooth(flux,3),col=220
 wait,1                            
 endfor 
goto,start
end

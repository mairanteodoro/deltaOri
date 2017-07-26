function tvs, list_sp, list_sn, fout1, sig0=sigma0, data

!except=0

;+ TVS.PRO
; Baseado no paper de Fullerton et al. 1996
;
; MODO DE USAR:
;  - Interativo:
;    	result = tvs()
;  - Direto:
;   	result = tvs('list', 'sn', 'fileout')
;   	onde
;   	'list' eh a lista dos espectros;
;   	'sn' eh a lista com os respectivos S/N;
;   	'fileout' eh o arquivo que vai conter
;   	          os resultados;
;   	'sigma0' eh a variavel que vai conter
;   	    	  o valor do desvio padrao da
;   	    	  distribuicao do chi quadrado.
;
; OUTPUT:
;   result -> matriz [npix, 3] contendo
;   	      lambda 	      -> result[*,0];
;   	      tvs   	      -> result[*,1];
;   	      espectro medio. -> result[*,2];
;   Alem disso, o resultado tambem
;   eh salvo no arquivo denominado
;   'fileout'.
;
;  N.B.: a distribuicao do (TVS)^1/2 segue a
;        distribuicao do chi^2, isto eh, o erro
;        a ser utilizado como significancia
;        estatistica eh o seguinte:
;
;        (TVS)^1/2 ~ (chi^2)^1/2 x sigma0
;
;        onde sigma0 eh o valor retornado por
;        esta tarefa. Para obter as curvas de
;        significancia estatistica, basta 
;        substituir o valor de chi^2 por
;        1 (68.5%), 2 (95%) ou 3 (99.9%).
;
; CRIADO EM 21/11/2008
;
; MODIFICACOES:
;

on_error, 2
if n_params() lt 3 then begin

    list_sp = '' ; lista dos espectros, um por linha
    list_sn = '' ; lista contendo os respectivos S/N de cada espectro
    fout1 = '' ; arquivo contendo duas colunas: lambda e TVS
    print, '++++++++++++++++++++++++++++++++++++++++++++'
    print, ''
    print, 'O numero de parametros fornecidos eh < 3.'
    print, '     Entrando no modo interativo'
    print, ''
    read, list_sp, prompt='- Lista com os espectros: '
    read, list_sn, prompt='- Lista com o S/N de cada espectro: '
    print, ''
    read, fout1, prompt='- Nome do arquivo de saida para os resultados: '
    print, ''

endif

print, ' Lendo os arquivos de entrada:'
print, '   -> ', strcompress(list_sp, /rem)
print, '   -> ', strcompress(list_sn, /rem)

files = rd_tfile(list_sp, /quiet)
sn    = rd_tfile(list_sn, /conv, /quiet)

; numero de espectros (n) e pixels (npix)
n    = n_elements(files)
npix = file_lines(files[0])

print, 'n/npix = ', strcompress(n, /rem), '/', strcompress(npix, /rem)

; definicao:
; i -> espectro
; j -> lambda do i-esimo espectro
;**************

; extracao das colunas lambda e intensidade
; para facilitar a leitura de um ou de outro
; separadamente:
openw, lun_lambda, 'lambda.txt', /get_lun
openw, lun_intens, 'intens.txt', /get_lun
for i=0, n-1 do begin
    temp = rd_tfile(files[i], 2, /conv, /quiet)
    printf, lun_lambda, temp[0,*]
    printf, lun_intens, temp[1,*]
endfor
close, /all

; criacao da matriz S:
print, ' Calculando a matriz S...'
s = make_array(n, npix, /double)
s_all = rd_tfile('intens.txt', /conv, /quiet)
k = 0 ; contador
for i=0, n-1 do begin
    for j=0, npix-1 do begin
    	s[i,j] = s_all[k]
		k++
    endfor
endfor

; guardar o valor de lambda:
l = make_array(n, npix, /double)
l_all = rd_tfile('lambda.txt', /conv, /quiet)
for i=0, n-1 do begin
    for j=0, npix-1 do begin
    	l[i,j] = l_all[j]
    endfor
endfor

; [..] for strong exposures,
; read-out noise can be neglected
; and
alfa = s

sn0 = sqrt(total(sn^2)/n); S/N medio (eq. 10)
snc = sn ; S/N no continuo de cada espectro

; dispersao do continuo em torno
; da sua media eh dado por:
sigma0 = 1./sn0 ; medio
sigmac = 1./snc ; continuo

; peso de cada espectro
; baseado no seu S/N:
w = (sigma0/sigmac)^2 ;(eq. 11)

; matriz S barra (eq. 14):
print, ' Calculando o espectro medio (Sbarra)...'
s_bar = make_array(npix, /double)
for j=0, npix-1 do begin
    s_bar[j] = total(w*s[*,j])
endfor
s_bar = s_bar/n

; matriz diferenca (eq. 13):
print, ' Calculando a matriz diferenca (D)...'
d = make_array(n, npix, /double)
for i=0, n-1 do begin
    for j=0, npix-1 do begin
    	d[i,j] = sqrt(w[i]/alfa[i,j])*(s[i,j]-s_bar[j])
    endfor
endfor

; temporal variance spectrum (eq. 4):
print, ' Calculando o TVS...'
tvs = make_array(npix, /double)
for j=0, npix-1 do begin
    tvs[j] = total(d[*,j]^2)
endfor
tvs = tvs/(n-1.)

print, ''
print, '-> Salvando os resultados em ', strcompress(fout1, /rem), '...'
openw, lun_tvs, fout1, /get_lun
for j=0, npix-1 do begin
    if j eq 0 then begin
    	printf, lun_tvs, '####################################################'
	printf, lun_tvs, '# ', systime()
    	printf, lun_tvs, "# A lista de entrada '", strcompress(list_sp,/rem), "'"
	printf, lun_tvs, '# continha os seguintes espectros:'
	printf, lun_tvs, '#'
	for i=0, n-1 do begin
    	    printf, lun_tvs, '# ', strcompress(i+1,/rem),$
	    	' - ', strcompress(files[i],/rem),$
		' (S/N=',  strcompress(sn[i],/rem),')'
	endfor
	printf, lun_tvs, '#'
    	printf, lun_tvs, '# S/N(medio)=', strcompress(sn0, /rem)
    	printf, lun_tvs, '# sigma(medio)=', strcompress(sigma0, /rem)
	printf, lun_tvs, '#'
	printf, lun_tvs, '#+++++++++++++++++++++++++++++++++++++++++++++++++++'
	printf, lun_tvs, "# N.B.: caso existam '-NaN', utilizar o seguinte"
	printf, lun_tvs, '#       para plotar:'
	printf, lun_tvs, "#-> tvs = rd_tfile('tvs.txt', 3, /conv)"
	printf, lun_tvs, '#-> x  = tvs[0, where(finite(tvs[1,*], /nan) ne 1)]'
	printf, lun_tvs, '#-> y1 = tvs[1, where(finite(tvs[1,*], /nan) ne 1)]'
	printf, lun_tvs, '#-> y2 = tvs[2, where(finite(tvs[1,*], /nan) ne 1)]'
	printf, lun_tvs, '#+++++++++++++++++++++++++++++++++++++++++++++++++++'
	printf, lun_tvs, '#'
	printf, lun_tvs, '# (lambda) (TVS) (espectro_medio)'
    	printf, lun_tvs, '####################################################'
    endif
    printf, lun_tvs, strcompress(l[0,j], /rem), ' ', $
    	    	     strcompress(tvs[j], /rem), ' ', $
		     strcompress(s_bar[j], /rem)
endfor
close, /all

; limpeza:
file_delete, ['lambda.txt','intens.txt'], /allow, /quiet

print, ''
print, '+++++++++++++++++FEITO++++++++++++++++++++++'

data = [ [reform(l[0,*])], [tvs], [s_bar] ]

return, data


end

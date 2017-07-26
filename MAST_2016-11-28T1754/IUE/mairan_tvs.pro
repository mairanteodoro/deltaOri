pro mairan_tvs

; ####################################################
; PROPOSITO:
; Programa para calculo do temporal variance spectrum.
;
; MODO DE USAR:
; mairan_tvs
;
; ENTRADAS:
; -> lista dos espectros a serem analisados
;    ('dispcorrected' e normalizados)
; -> lista do sinal-ruido de cada espectro
;
; SAIDA:
; Os espectros sao comparados pixel a pixel
;  e o espectro TVS eh salvo em um arquivo
;  no formato [lambda, tvs].
;
; DATA DA CRIACAO DESTE PROG:
; 25/06/2008
;
; MODIFICACOES:
;
; 18/11/2008:
; - O programa agora coloca o valor do S/N medio
;   na primeira linha do arquivo de TVS.
;
; ####################################################

on_error, 2
!except=0

; INPUTS:
dirin='' & list='' & sn = ''
read, prompt='DIR (in/output): ', dirin
read, prompt='LIST (one spec per line): ', list
read, prompt='S/N (one per line): ', sn

; OUTPUTS:
dirout = dirin
fileout1 = 'lambda'
fileout2 = 'intens'


;###################
; PARTE 1 (tvs_1.c)
;###################
nlines = file_lines(dirin+list)
spec_list = rd_tfile(dirin+list)

openw, lun1, dirout+fileout1, /get_lun
openw, lun2, dirout+fileout2, /get_lun
print, 'LENDO E PREPARANDO OS ARQUIVOS...'
for i=0, nlines-1 do begin
    data = rd_tfile(dirin+spec_list[i], 2, /convert)

    ; escrever os lambdas
    printf, lun1, data[0,*]

    ; escrever as intensidades
    printf, lun2, data[1,*]
endfor
close, lun1 & free_lun, lun1
close, lun2 & free_lun, lun2
print, ''
print, 'FEITO!'
print, ''


;###################
; PARTE 2 (tvs_2.c)
;###################

print, 'CALCULANDO...'
n = double(nlines) ; numero de espectros na lista
npix = size(data)
pix = double(npix[2]) ; numero de pixels do espectro

sn = double(rd_tfile(dirin+sn, /convert)) ; S/N dos espectros

; calcular o S/N medio
snmedio = 0.
for i=0, n-1 do begin
    snmedio = snmedio + sn[i]^2
endfor
snmedio = snmedio / n
snmedio = sqrt(snmedio)

; leitura das intensidades dos espectros
intens = double(rd_tfile(dirin+fileout2))
; numero de pontos do espectro
nspec = file_lines(dirin+spec_list[0])

med = make_array(pix, /double)
fileout_med = 'intens_average.txt'
openw, lun_med, dirout+fileout_med, /get_lun
print, '>>>>>>>> INTENSIDADES MEDIAS...'
; EQUACAO (14) DO FULLERTON ET AL. 1996
for p=0., pix-1. do begin
    for q=0., n-1. do begin
    	value = intens[q*pix+p]
	med[p] = med[p] + value * (sn[q] / snmedio)^2.
    endfor
    med[p] = med[p] / n
    printf, lun_med, med[p]
endfor
close, lun_med & free_lun, lun_med



; leitura das intensidades dos espectros
lambda = rd_tfile(dirin+fileout1)

lmed = make_array(pix, /double)
fileout_lam = 'lamb_average.txt'
openw, lun_lam, dirout+fileout_lam, /get_lun
print, '>>>>>>>> COMPRIMENTOS DE ONDA MEDIOS...'
for p=0., pix-1. do begin
    for q=0., n-1. do begin
    	value3 = lambda[q*pix+p]
	lmed[p] = lmed[p] + value3
    endfor
    lmed[p] = lmed[p] / n
    printf, lun_lam, lmed[p]
endfor
close, lun_lam & free_lun, lun_lam


diff = make_array(n,pix,/double)
; calcula a matriz diferenca
fileout_diff = 'diferencas.txt'
openw, lun_diff, dirout+fileout_diff, /get_lun
print, '>>>>>>>> MATRIZ DIFERENCA...'
; EQUACAO (13) DO FULLERTON ET AL. 1996
for s=0., n-1. do begin
    for t=0., pix-1. do begin
    	diff[s,t] = (sqrt(sn[s]/snmedio))*(intens[s*pix+t] - med[t]);/(sqrt(intens[s*pix+t]))
	printf, lun_diff, diff[s,t]
    endfor
endfor
close, lun_diff & free_lun, lun_diff


tvs = make_array(pix, /double)
; calcula o TVS *****************************
fileout_tvs = 'tvs.txt'
openw, lun_tvs, dirout+fileout_tvs, /get_lun
print, '>>>>>>>> TEMPORAL VARIANCE SPECTRUM...'
; EQUACAO (4) DO FULLERTON ET AL. 1996
for p=0, pix-1 do begin
    tvs[p] = 0
    for q=0, n-1 do begin
    	value2 = diff[q,p]
	tvs[p] = tvs[p] + value2^2
    endfor
    tvs[p] = tvs[p] / (n-1)
    if p eq 0 then begin
    	printf, lun_tvs, '# S/N medio: '+strcompress(string(snmedio), /remove)
    	printf, lun_tvs, lmed[p], med[p], tvs[p]
    endif
    printf, lun_tvs, lmed[p], med[p], tvs[p]
endfor
close, lun_tvs & free_lun, lun_tvs
print, 'FEITO!'

print, ''
print, '-> ARQUIVO CONTENDO O TVS: '+strcompress(dirout+fileout_tvs,/remove)

files = ['intens_average.txt','lamb_average.txt','diferencas.txt']
file_delete, files, /ALLOW_NONEX, /QUIET
print, ''



end


************************************
******MODELO LINEAL GENERAL*********
*===================================
* Heber Baldeon
* Grupo Lambda
*===================================
************************************

/*The Statistical Software Components (SSC) archive, 
hosted by Boston College (USA) and maintained by Christopher F. Baum, 
is the largest collection of community-contributed Stata programs 
for data manipulation, statistics, and graphics. 
Existing Stata commands search, ssc describe, and ssc install 
make locating and installing community-contributed programs a snap.*/

ssc whatshot
ssc whatshot, n(15)
ssc whatshot, author(cox) 
ssc whatshot, author(baum) 

*I. Subir la información

*cd
cd "C:\Users\percy\OneDrive\Documentos\ECONOMETRIA APLICADA - LAMBDA\Sesion 01"
use Ec_Ingreso.dta

*Corregir las etiquetas de las variables
label variable pet "Población en edad de trabajar"
label variable edu007 "Número de años de educación aprobados"
label variable exp "Años de experiencia"
label variable exp2 "Años de experiencia al cuadrado"
label variable fex "Factor de Expansión"
label variable mig "Inmigrantes en los últimos 5 años"
label variable union "Tiene algun tipo de unión conyugal"
label variable men5 "Número de menores de 5 años en el hogar"

label define sexo 1 "Hombre" 0 "Mujer"
label values e03 sexo	


*II. Analisis Preliminar (descriptivo)

browse // Para ver la data
describe 
summarize // Resumen estadístico descriptivos
codebook e03 // Resumen estadístico y de la data (missing values)
codebook lny edu007 

*d, su, codebook

sum edu007 e03
table edu007 e03 
tab edu007 e03 
summ lny lnm edu007 exp exp2 e03, detail

*Correlación de Pearson (relación lineal)
cor lny lnm edu007 exp exp2
*Covarianza
cor lny lnm edu007 exp exp2, cov

*Correlación pairwise (dependencia en segundo momento o cuadrática)
pwcorr lny lnm edu007 exp exp2
pwcorr lny lnm edu007 exp exp2, sig
pwcorr lny lnm edu007 exp exp2, star(.01) bonferroni obs // Bonferroni: Ajustes estadísticos para evaluar la significancia

*Correlación de Spearman y Tau-Kendall (Test no paramétricos para diferentes tipos de relación)
spearman lny lnm edu007 exp exp2
ktau lny lnm edu007 exp exp2

graph matrix lny lnm edu007 exp exp2, name(G1, replace)
graph matrix lny exp exp2, by(e03) name(G2, replace)

twoway (scatter lny lnm, sort)
twoway (scatter lny edu007, sort)
twoway (scatter lny exp, sort)
twoway (scatter lny exp2, sort)

egen lny_ed=mean(lny), by(edu007) 
label variable lny_ed "Ingreso Promedio por año de educacion"

*Relación entre los ingresos y el nivel de estudio
sum lny_ed
*Summarize guarda los escalares con r: r(N), r(mean), r(p50)
global tlny_ed=r(mean)
*local
global var_3 lnm edu007 exp exp2
display $tlny_ed

*reg lny $var_3

sum edu007
global edu0071=r(mean)
display $tedu007

line lny_ed edu007
line lny_ed edu007, sort name(G3, replace)
scatter lny_ed edu007, sort xline($tedu007) yline($tlny) name(G4, replace)
 
*III. Estimación

regress lny lnm
predict lny_1, xb
twoway (scatter lny lnm, sort) (line lny_1 lnm, sort)

summ lny lnm edu007 exp exp2

regress lny lnm edu007 exp exp2  
est store MCO_1
estat ic

*Estimadores MELI (el MCO es más eficiente que el GMM porque es MELI)

gmm (lny-lnm*{b1}-edu007*{b2}-exp*{b3}-exp2*{b4}-{b5}), instruments(lnm edu007 exp exp2) nolog
gmm (lny-{xb:lnm edu007 exp exp2}-{b0}), instruments(lnm edu007 exp exp2) nolog
est store GMM_1

estimate table MCO_1 GMM_1, b se

*Pruebas de Hipótesis lineal

test lnm=1
test (lnm) (edu007+exp=1)
test lnm=edu007=exp=exp2=0

*Pruebas de hipótesis no lineal (cambia la naturaleza del estadístico)
testnl _b[lnm]/_b[edu007]=1
test lnm=edu007

*Cambio de pendiente

xi: regress lny lnm exp exp2 i.edu004
test _Iedu004_2 _Iedu004_3 _Iedu004_4 _Iedu004_5 _Iedu004_6
estat sum

xi: regress lny lnm exp exp2 i.edu004*edu007  
testparm _Iedu004_2 - _Iedu004_6 _IeduXedu00_2 - _IeduXedu00_6

*Test de cambio de régimen (Solo aplicable a series de tiempo)
*Cuando conoces el periodo del cambio de régimen
*estat sbknown, break(t)

*IV. Evaluación de los supuestos del MCO

*1. Multicolinealidad

graph matrix lnm edu007 exp exp2, half name(G5, replace)

regress lny lnm edu007 exp exp2  
vif

*2. Homocedasticidad

*Metodo Gráfico (residual vs fitted)
rvfplot, yline(0) name(G6, replace)
rvfplot, yline(1) xline(10) name(G7, replace)

regress lny lnm edu007 exp exp2
predict plny

gen error1=lny-plny
predict error2, resid

scatter error1 plny, yline(0) name(G8, replace)

*Metodo Formal

hettest 
imtest, white

*Si solo hay problemas de heteroscedasticidad
regress lny lnm edu007 exp exp2, robust   //{n/(n-k)}(u_j)^2
*Alternativamente
regress lny lnm edu007 exp exp2, vce(robust)

*Si hay problemas de heteroscedasticidad y correlación entre errores
regress lny lnm edu007 exp exp2, vce(cluster edu007)
estimate store MCOrobust

*Otras alternativas
regress lny lnm edu007 exp exp2, vce(bootstrap)
estimate store MCOboots
regress lny lnm edu007 exp exp2, vce(hc2)  //{n/(n-k)}(u_j)^2/(1-h_jj)
estimate store MCOhc2
regress lny lnm edu007 exp exp2, vce(hc3)  //{n/(n-k)}(u_j)^2/(1-h_jj)^2
estimate store MCOhc3

estimate table MCOrobust MCOboots MCOhc2 MCOhc3, se stats(N r2 F ll)

*Mínimos cuadrados generalizados

vwls lny lnm edu007 exp exp2, sd(edu007)
glm lny lnm edu007 exp exp2, family(gaussian) link(identity)

*3. Autocorrelacion (aplicable a series de tiempo)

*tsset 
estat dwatson
prais vardep varindep, corc
newey vardep varindep, lag(2)

*4. Normalidad de los residuos

predict res_mod1, r 
predict lny_est, xb 

*Metodo Gráfico 

histogram res_mod1, name(G9, replace)

qnorm res_mod1, name(G10, replace)

kdensity res_mod1, name(G11, replace)
kdensity res_mod1, normal name(G12, replace)

*Metodo Formal (no aplicable con opcion robust o modelo corregido)

regress lny lnm edu007 exp exp2
  
predict lnyr, residual
predict lnyrt, rstandard
predict lnyrs, rstudent

summarize lnyr*

sktest lnyr-lnyrs
swilk lnyr-lnyrs
sfrancia lnyr-lnyrs

*5. Variables redundantes (Sesgo de especificación)

*Ojo revisar si es sesgo de especificación o inconsistencia

regress lny edu007 exp exp2  
est store reg1

regress lny lnm edu007 exp exp2
est store reg2

hausman reg1 reg2

*6. Variables omitidas (Sesgo de especificación)

regress lny lnm edu007 exp exp2 e03
estat ovtest 

*7. Bondad de ajuste

predict lnyp if e(sample), xb

egen plny_ed=mean(lnyp), by(edu007) 

line plny_ed lny_ed edu007, sort  name(G13, replace)

*8. Estimacion con el factor de expansion

regress lny lnm edu007 exp exp2
est store a1
regress lny lnm edu007 exp exp2 [iw=fex]
est store a2
regress lny lnm edu007 exp exp2 [pw=fex], robust
est store a3

estimate table a1 a2 a3
estimate table a1 a2 a3, se stats(N r2 F ll) b(%10.4f) title(Comparación de modelos)

*9. Presentación

*Instalación
ssc install outreg2

reg lny lnm edu007 exp exp2 if e03==1 
outreg2 using primero.xls, ctitle(Hombre) 
reg lny lnm edu007 exp exp2 if e03==0 
outreg2 using primero.xls, append ctitle(Mujer) 

*latex

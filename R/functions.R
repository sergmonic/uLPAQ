
# Función interna para factores de normalización de distribuciones comunes
factoresEstadarizacion<-function(distrib){
  distribuciones=data.frame("rectangular"=sqrt(3),"triangular"=sqrt(6))
  distribuciones[[distrib]]
}

#' Incertidumbre estándar para material volumétrico de vidrio y micropipetas tipo A
#'
#' Devuelve la incertidumbre asociada a intrumentos volumétricos (pipetas, buretas, balores aforados o micropipetas) empleado en el LPAQ, teniendo
#' en cuenta las siguientes fuentes principales:
#'  1. Error máximo permitido u(V,EM), considerando una distribucion triangular o calibración (sesgo, precisión) u(V,sesgo), u(V,R)
#'  2. Desviación de temperatura respecto a 20 ºC, u(V,T)
#'
#' Los valores devueltos y las especificaciones empleadas se expresan en mL, teniendo en cuenta los valores de las normas
#' ISO 385:2005 (Buretas) (1), ISO 648:2008 (Pipetas volumétricas) (2), ISO 1042:1998 (Balones volumétricos) (3) y
#' ISO 8655-2:2002 (micropipetas) (4).
#'
#' Si el valor de incertidumbre del sesgo
#'
#' @param nominal valor nominal del material
#' @param tipo tipo de material volumétrico: "pipeta", "bureta", "balon", "micropipeta"
#' @param clase A o AS: "A_AS"), B: "B" (Solo para material de vidrio)
#' @param temp temperatura de uso del instrumento volumétrico
#' @param subdivision solo para buretas 0.01 mL, 0.02 mL, 0.05 mL, 0.1 mL y 0.2 mL
#' @param volumen entregado para el caso de buretas y la micropipetas. Los valores se deben expresar en "mL"
#' @param calibracion vector  c(sesgo, u_s, prec,k_us), con sesgo (error); u_s, incertidumbre estándar del sesgo, incertidumbre estándar debida a la precisión u(V,Prec) y
#' factor de cobertura de la incertidumbre del sesgo (k_us=1.0 por omisión). Si el valor de la especificación corresponde a "Error %", "Inaccuracy %" o "Accuracy %"
#' se recomienda usar k_us=sqrt(3), asumiendo un distribución rectangular para este parámetro.
#' La incertidumbre u(V,Prec) se relaciona con la variación del volumen en el llenado del instrumento (depende del líquido y del operador), esta se puede expresar como una desviación estándar en el caso de una única medición.
#' Los valores se deben expresar en "mL"
#' @param gamma, coefficiente de expansión t\'ermico (por omisión gamma=2.1E-4/ºC para el agua)
#' @param distrib distribución asumida (por omisión distrib="triangular")
#' @param t_ref Temperatura de referencia a la que se caracterizó el instrumento (20 ºC por omisión)
#'
#' @return data.frame con volumen nominal o corregido, incertidumbre estandar combinada (u_V), incertidumbre
#' debida a la calibración o al error máximo permitido (u_V_c) e incertimubre debida a variación de la temperatura (u_V_T).
#' @examples
#' instVolumetrico(nominal=100,tipo="balon",clase="B")
#'
#' @references
#' (1) ISO 385:2005, Laboratory glassware - Burettes. Geneva, Switzerland. \url{https://www.iso.org/standard/38678.html}
#'
#' (2) ISO 648:2008, Laboratory glassware - Single-volume pipettes. Geneva, Switzerland. \url{https://www.iso.org/standard/44142.html}
#'
#' (3) ISO 1042:1998, Laboratory glassware - One-mark volumetric flasks. Geneva, Switzerland. \url{https://www.iso.org/standard/25484.html}
#'
#' (4) ISO 8655-2:2002, Piston-operated volumetric apparatus — Part 2: Piston pipettes. Geneva, Switzerland. \url{https://www.iso.org/standard/29727.html}
#'
#' @export
instVolumetrico<-function(nominal,tipo,clase="A_AS",subdivision=0.1,volumen=0,temp=20.0,calibracion=c(0.0,0.0,0.0,1.0),gamma=2.1E-4,distrib="triangular",t_ref=20.0){
  .tipo=c("bureta"="buretas","pipeta"="pipetas","balon"="balones","micropipeta"="micropipetas")[[tipo]]
  k_us=calibracion[4]
  V=nominal
  if(.tipo=="buretas"){
    if(volumen>0.0){V=volumen}
    EMP=EMP_materialVolumetrico[["buretas"]][EMP_materialVolumetrico[["buretas"]]$Nominal==nominal &
                                     EMP_materialVolumetrico[["buretas"]]$Subdivision==subdivision,][[clase]]
  }else{
    if(.tipo=="micropipetas"){
      if(volumen>0.0){V=volumen}
      EMP=EMP_materialVolumetrico[["micropipetas"]][EMP_materialVolumetrico[["micropipetas"]]$Nominal==nominal*1000,]
      EMP_s=EMP[[3]]/1000
      EMP=EMP[[2]]/1000
    }else{
      EMP=EMP_materialVolumetrico[[.tipo]][EMP_materialVolumetrico[[.tipo]]$Nominal==nominal,][[clase]]
    }
  }

  V=V-calibracion[1]

  # u debida a veracidad o a datos de calibracion
  if(calibracion[2]<=1E-15){
    if(calibracion[3]>=1E-15){
      u_V=sqrt((EMP/factoresEstadarizacion(distrib))^2+calibracion[3]^2)
    }else{
      if(.tipo=="micropipetas"){
        u_V=sqrt((EMP/factoresEstadarizacion(distrib))^2+EMP_s^2)
      }else{
        u_V=EMP/factoresEstadarizacion(distrib)
      }
    }
  }else{
    u_V= sqrt((calibracion[2]/k_us)^2+calibracion[3]^2)
  }

#Efecto de la temperatura
u_VeT=u_VT(V,DT = abs(temp-t_ref),gamma = gamma)
u_Vt=sqrt(u_V^2+u_VeT^2)

data.frame(V=V,u_V=u_Vt,u_V_c=u_V,u_V_T=u_VeT,EMP=EMP)

}


#' Incertidumbre en el volumen de un líquido debido a variaciones de temperatura
#'
#' Devuelve la incertidumbre en el volumen del líquido empleado, debido a variaciones de temperatura  respecto a la temperatura estándar de referencia de 20 ºC.
#' No se considera la expansión del material del instrumento, dado que los coeficientes de expansión térmica del vidrio son más bajos que los del líquido (Borosilicato 3.3 9.9E-6/ºC
#' y vidrio común 27E-6/ºC). El efecto de la temperatura sobre el vidrio produce variaciones en volumen de 0.007 % para el borosilicato y 0.02 % para el vidrio común
#' a una temperatura de 27 ºC, evidenciando un  efecto es mínimo de la temperatura sobre la expansión del material del instrumento (1).
#'
#' @param volumen Volumen utilizado o entregado en la medición
#' @param DT temperatura respecto a la temperatura estándar de referencia (normalmente 20 ºC) en grados Celsius
#' @param gamma, coefficiente de expansión térmico (por omisión gamma=2.1E-4/ºC para el agua)
#' @param distrib distribución asumida (por omisión distrib="rectangular")
#'
#' @return incertidumbre estándar en el volumen del líquido debido a variación de temperatura respecto a la referencia.
#'         La incertidumbre calculada se obtiene asumiendo una distribución rectangular para la variación de temperatura.
#' @examples
#' u_VT(volumen=100.0,DT=4)
#'
#' @references
#' (1) ISO 4787:2010, Laboratory glassware - Volumetric instruments - Methods for testing of capacity and for use. Geneva, Switzerland. \url{https://www.iso.org/standard/41807.html}
#'
#' @export
u_VT<-function(volumen,DT, gamma=2.1E-4,distrib="rectangular"){
  volumen*gamma*DT/factoresEstadarizacion(distrib)
}


#' Incertidumbre en mediciones de masa, empleando especificaciones del fabricante (balanzas sin calibración)
#'
#' Devuelve la incertidumbre debida a una operación de pesaje en una balanza de la que no se tiene información de
#' calibración. Teniendo en cuenta las siguientes consideraciones:
#'    i)   La función solo considera las balanzas que se emplean en el LPAQ (ME204, AL54 y AS220)
#'    ii)  La muestra y el recipiente en que se mide, se colococan centrados en el plato (i.e. no se considera el efecto de la excentricidad)
#'    iii) La muestra tiene masa definida i.e. efectos de volatilidad de la muestra o higroscopicidad son despreciables
#'    iv)  El operador de la balanza conoce el procedimiento de pesaje y lo efectua con el debido cuidado
#'
#'
#' @param lectura Lectura de masa en unidad "g", obtenida del instrumento de medición.
#' @param densidad densidad de la muestra en kg/m^3
#' @param balanza empleada:"ME204" (valor por omisión),"AL54" y AS220
#' @param u_densidad incertidumbre estándar relativa en la densidad de la muestra, 5 % por omisión y tomada como una distribución rectangular.
#' Si se conoce la incertidumbre estándar de la densidad, se debe incluir como u_den*sqrt(3)
#' @param T_ c(T_prom,T_min,T_max) temperatura del aire en ºC
#' @param p c(p_prom,p_min,p_max) presión expresada en hPa
#' @param HR c(HR_prom,HR_min,HR_max) porcentaje de humedad relativa
#' @param den_cal c(denCal,u_denCal) densidad de las pesas empleadas para la calibración e incertidumbre asociada en kg/m^3,
#' por omision den_cal=c(8006,10)
#' @param den_aire_cal c(denAireCal,u_denAireCal) densidad e incertidumbre estándar del aire durante el ajuste de la sensibilidad de la balanza en kg/m^3,
#' por omision den_aire_cal=c(1.2,0.0129)
#' @param x_CO2 fracción molar de CO2 (por omisión 0.00040)
#'
#' @return lista con:
#'   i)   data.frame (ms) para masa corregida por efecto de flotación (ms) e incertidumbre estándar en la medición de masa (u_ms) expresada en "g"
#'   ii)  data.frame (Bu) con corrección por efecto de flotación (adimensional), y factores de influencia considerados: a) densidad de la muestra (d_s), b) densidad del aire durante la medición (d_as), c) densidad del aire durante la ajuste y calibración (d_ac), c) densidad de las pesas empleandas para ajuste y calibración (d_pc)
#'   iii) data.frame (Ws) con incertidumbre debidad al instrumento  (ws) en "mg", y sus factores de influencia: a) Resolución (u_Res), b) repetibilidad (u_rep), c) desvio de la linealidad (u_NL), d) Tolerancia en la sensibilidad (u_TS), e) Efecto térmico en la sensibilidad (u_ETS)
#' @examples
#' u_masa(.0010,1636)
#'
#' @references
#' (1) Reichmuth, A., Wunderli, A., Weber, M., Meyer, V.R., The Uncertainty of Weighing Data Obtained with Electronic Analytical Balances,
#' Microchim. Acta 148, 133–141 (2004). \url{https://doi.org/10.1007/s00604-004-0278-3}
#'
#' @export
u_masa<-function(lectura,densidad,balanza="ME204",u_densidad=0.05,T_=c(21.5,15.7,24.7),p=c(749.6,746.0,753.8),HR=c(47.0,31.2,67.4),den_cal=c(8006,10),den_aire_cal=c(1.2,0.0129),x_CO2=0.00040){

m_s_mg=lectura*1000 # Expresa la lectura en mg
dT=(T_[3]-T_[2])/2

#=========================================================================
#Efecto de flotación por empuje del aire
#=========================================================================
d_s=c(densidad,u_densidad*densidad/sqrt(3))
Bu_s=flotacion(den_s = d_s,T_ = T_,p=p,HR = HR,den_cal = den_cal,den_aire_cal = den_aire_cal,x_CO2 = x_CO2)


#=========================================================================
#Efectos propios del instrumento: cada componente se expresada en mg
#=========================================================================

# 1. Resolución de la balanza, considerando la operación de ajuste de cero o tara, de acuerdo con la forma de pesado
u_Res=balanzasLPAQ[[balanza]][2]/sqrt(6) #mg

# 2. Repetibilidad
u_rep=balanzasLPAQ[[balanza]][3] # mg

# 3. Desviación de la linealidad
u_desvLinealidad=2.0*balanzasLPAQ[[balanza]][4]/sqrt(3) # mg

# 4. Tolerancia en la pendiente lectura/masa (Sensibility Tolerance)
u_ST=lectura/1000*balanzasLPAQ[[balanza]][6]/sqrt(3) # mg

# 5. Variación térmica de la sensibilidad
u_TC=lectura/1000*balanzasLPAQ[[balanza]][5]/sqrt(3)*dT # mg

u_lectura=sqrt(u_Res^2+u_rep^2+u_desvLinealidad^2+u_ST^2+u_TC^2)

######COMBINACION DE FUENTES DE INCERTIDUMBRE

ms_=Bu_s$Bu*lectura
u_ms=sqrt((Bu_s$u_Bu/Bu_s$Bu)^2+(u_lectura/m_s_mg)^2)*ms_

ms=data.frame(ms=ms_,u_ms=u_ms)
bu=data.frame(Bu=c(Bu_s$Bu,Bu_s$u_Bu),d_s=d_s,d_as=c(Bu_s$das,Bu_s$u_das),d_ac=den_aire_cal,d_pc=den_cal)
ws=data.frame(ws=m_s_mg,u_ws=u_lectura,u_Res,u_rep,u_NL=u_desvLinealidad,u_TS=u_ST,u_ETS=u_TC)

list("m_s"=ms,"Bu"=bu,"Ws"=ws)

}


#' Minima cantidad a pesar
#'
#' Devuelve mínima cantidad a pesar para una incertidumbre estandar de pesada requerida
#'
#' @param u_r  incertidumbre estándar relativa requerida, expresada porcentualmente.
#' @param densidad densidad de la muestra en kg/m^3
#' @param balanza empleada:"ME204" (valor por omisión),"AL54" y AS220
#' @param u_densidad incertidumbre estándar relativa en la densidad de la muestra, 5 % por omisión y tomada como una distribución rectangular.
#' Si se conoce la incertidumbre estándar de la densidad, se debe incluir como u_den*sqrt(3)
#' @param T_ c(T_prom,T_min,T_max) temperatura del aire en ºC
#' @param p c(p_prom,p_min,p_max) presión expresada en hPa
#' @param HR c(HR_prom,HR_min,HR_max) porcentaje de humedad relativa
#' @param den_cal c(denCal,u_denCal) densidad de las pesas empleadas para la calibración e incertidumbre asociada en kg/m^3,
#' por omision den_cal=c(8006,10)
#' @param den_aire_cal c(denAireCal,u_denAireCal) densidad e incertidumbre estándar del aire durante el ajuste de la sensibilidad de la balanza en kg/m^3,
#' por omision den_aire_cal=c(1.2,0.0129)
#' @param x_CO2 fracción molar de CO2 (por omisión 0.00040)
#' @param polyOrden grado de polinomio usado para ajustar la curva caracterítica de la balanza
#'
#' @return Masa minima a pesar para la incertidumbre estándar especificada y curva caracterítica de la balanza
#'
#' @importFrom stats lm
#' @importFrom latex2exp TeX
#'
#' @references
#' (1) Reichmuth, A., Wunderli, A., Weber, M., Meyer, V.R., The Uncertainty of Weighing Data Obtained with Electronic Analytical Balances,
#' Microchim. Acta 148, 133–141 (2004). \url{https://doi.org/10.1007/s00604-004-0278-3}
#'
#' @export
masa_minima<-function(u_r,densidad,balanza="ME204",u_densidad=0.05,T_=c(21.5,15.7,24.7),p=c(749.6,746.0,753.8),HR=c(47.0,31.2,67.4),den_cal=c(8006,10),den_aire_cal=c(1.2,0.0129),x_CO2=0.00040,polyOrden=4){
  # u_r=0.1
  # densidad=1634
  # balanza="XPE205"
  # u_densidad=0.05
  # T_=c(21.5,15.7,24.7)
  # p=c(749.6,746.0,753.8)
  # HR=c(47.0,31.2,67.4)
  # den_cal=c(8006,10)
  # den_aire_cal=c(1.2,0.0129)
  # x_CO2=0.00040
  # polyOrden=4
  #

  max_carga=balanzasLPAQ[[balanza]][1]*1000
  leg=balanzasLPAQ[[balanza]][2]
  x=seq(leg,max_carga,10)
  u_x=u_masa(x/1000,densidad=densidad,balanza=balanza,u_densidad=u_densidad,T_=T_,p=p,HR=HR,den_cal=den_cal,den_aire_cal=den_aire_cal,x_CO2=x_CO2)
  u_xr=u_x$m_s$u_ms/u_x$m_s$ms*100
  ylog=log(u_xr,10)
  xlog=log(x,10)
  model.pol=lm(ylog~poly(xlog,degree = polyOrden,raw = TRUE)) #Advertencia: Se ajusta con polinomio de grado 3, pero se debe evaluar en cada caso.
  coeffs.model=summary(model.pol)$coefficients[1:(polyOrden+1)]
  coeffs.model
  coeffs.model[1]=coeffs.model[1]-log(u_r,10)

  # plot(xlog,ylog-log(u_r,10))
  # xlog_=seq(min(xlog),max(xlog),0.01)
  # y=sapply(xlog_,function(x__,coeffs.model,order__){sum(coeffs.model*x__^seq(0,(order__)))},coeffs.model,order__)
  # lines(xlog_,y,col="blue")

  roots=polyroot(coeffs.model)
  roots=roots[which(Re(roots)>min(xlog) & Re(roots)<max(xlog) )]
  roots_indx=which(abs(Im(roots))<=1E-6)
  log_min_masa=Re(roots[roots_indx])
  min_mass=ceiling(10^log_min_masa[1])
  plot(x,u_xr,log="xy",type="b",pch=16,cex=0.5,col="blue",
       xlab="masa / mg", ylab=TeX(r'($u_{sr} / \%$)'),las=1, main=paste("Curva caracteristica de balanza ",balanza))
  abline(h=u_r,lty=2,col="red")
  abline(v=min_mass,lty=2,col="blue")
  text(min_mass,u_r,bquote(m[min]~"= "~.(min_mass)~"mg"),pos = 4)

  list("masa_min_mg"=min_mass,"curvaCaracteristica"=data.frame(masa_mg=x,u_sr=u_xr))

}



#' Densidad del aire
#'
#' Devuelve la densidad del aire a la condiciones ambientales especificadas, empleando la ecuación
#' revisada por el CIPM en 2007 (1). Esta ecuación es válida para 600 hPa <= p <= 1100 hPa y 15 ºC <= T <= 27 ºC.
#'
#' @param T_ temperatura del aire en ºC
#' @param p presión expresada en hPa
#' @param HR porcentaje de humedad relativa
#' @param x_CO2 fracción molar de CO2 (por omisión 0.00040)
#'
#' @return densidad del aire a las condiciones ambientales especificadas en kg/m^3
#' @examples
#' densidad_aire(T_=20,p=750,HR=50)
#'
#' @references
#' (1) Picard, A., Davis, R.S.,  Gläser, M., Fujii, K., Revised formula for the density of moist air (CIPM-2007), Metrologia, 45, 149-155 (2004) \url{https://doi.org/10.1088/0026-1394/45/2/004}.
#'
#' @export
densidad_aire<-function(T_,p,HR,x_CO2=0.00040){
 T_K=273.15+T_ # K
 p_Pa=p*100 # Pa
 h=HR/100 # fracción

 M_a=(28.96546+12.011*(x_CO2-0.0004))/1000 # kg/mol (Masa molar del aire seco)
 R_=8.314472 # J/(mol K) (Constante molar del gas ideal)
 M_v=18.01528/1000 # kg/mol (Masa molar del agua)

# Presión de vapor en el punto de saturación
#===================================
A_=1.2378847E-5 #K^-2
B_=-1.9121316E-2 # K^-1
C_=33.93711047
D_=-6.3431645E3 #K
p_sv=exp(A_*T_K+B_*T_K+C_+D_/T_K)

alpha=1.00062
beta=3.14E-8 # Pa^-1
gamma=5.6E-7 #K^-2

f_= alpha+beta*p_Pa+gamma*T^2


# Fracción molar de vapor de agua
#===================================

x_v=h*f_*p_sv/p_Pa

# Factor de compresibilidad
#============================
a0=1.58123E-6 #K Pa^-1
a1=-2.9331E-8 #Pa^-1
a2=1.1043E-10 # K^-1 Pa^-1
b0=5.707E-6 # K Pa^-1
b1=-2.051E-8 #Pa^-1
c0=1.9898E-4 #K Pa^-1
c1=-2.376E-6 #Pa^-1
d=1.83E-11 # K^2 Pa^-2
e=-0.765E-8 # K^2 Pa^-2

Z_=1-p_Pa/T_K*(a0 + a1*T_ + a2*T_^2 +(b0 + b1*T_)*x_v + (c0 +c1*T_)*x_v^2) + ((p_Pa/T_K)^2)*(d + e*x_v^2)

#densidad del aire
den_aire=p_Pa*M_a/(Z_*R_*T_K)*(1-x_v*(1-M_v/M_a))
den_aire
}

#' Efecto de flotación
#'
#' Devuelve el  factor de corrección debido al efecto de flotación (pesado en aire)
#' \deqn{Bu=\frac{1-\rho_{\mathrm{aire,cal}}/\rho_{\mathrm{pesas,cal}}}{ 1-\rho_{\mathrm{aire,med}}/\rho_{\mathrm{muestra}}}}
#'
#' @param den_s c(denMuestra,u_denMuestra) densidad e incertidumbre estándar de la muestra que se mide en kg/m^3
#' @param T_ c(T_prom,T_min,T_max) temperatura del aire en ºC
#' @param p c(p_prom,p_min,p_max) presión expresada en hPa
#' @param HR c(HR_prom,HR_min,HR_max) porcentaje de humedad relativa
#' @param den_cal c(denCal,u_denCal) densidad de las pesas empleadas para la calibración e incertidumbre asociada en kg/m^3,
#' por omision den_cal=c(8006,10)
#' @param den_aire_cal c(denAireCal,u_denAireCal) densidad e incertidumbre estándar del aire durante el ajuste de la sensibilidad de la balanza en kg/m^3,
#' por omision den_aire_cal=c(1.2,0.0129)
#' @param x_CO2 fracción molar de CO2 (por omisión 0.00040)
#'
#' @return data.frame con factor de corrección debido al efecto de flotación (Bu), incertidumbre estáncadar (u_Bu)
#' y densidad (das) e incertidumbre estándar del aire (u_das), en el ambiente del laboratorio.
#'
#' @examples
#' flotacion(c(998,30),p=c(1010,995,1025),T_=c(22,19,25),HR=c(50,25,75),den_aire_cal=c(1.19,0.0129))
#'
#' @references
#' (1) Picard, A., Davis, R.S.,  Gläser, M., Fujii, K., Revised formula for the density of moist air (CIPM-2007), Metrologia, 45, 149-155 (2004) \url{https://doi.org/10.1088/0026-1394/45/2/004}.
#'
#' @export
flotacion<-function(den_s,T_=c(21.5,15.7,24.7),p=c(749.6,746.0,753.8),HR=c(47.0,31.2,67.4),den_cal=c(8006,10),den_aire_cal=c(1.2,0.0129),x_CO2=0.00040){

den_aire_sup=densidad_aire(T_[2],p[3],HR[2])
den_aire_inf=densidad_aire(T_[3],p[2],HR[3])
den_aire_s=densidad_aire(T_[1],p[1],HR[1])
u_den_aire_s=(den_aire_sup-den_aire_inf)/sqrt(3)
Bu_=Bu(den_aire_s,den_aire_cal[1],den_cal[1],den_s[1])
ci_Bu=dBu(den_aire_s,den_aire_cal[1],den_cal[1],den_s[1])
u_Bu=sqrt(sum((ci_Bu*c(u_den_aire_s,den_aire_cal[2],den_cal[2],den_s[2]))^2))
data.frame(Bu=Bu_,u_Bu=u_Bu,das=den_aire_s,u_das=u_den_aire_s)
}

#=================================================================================
#Funcion auxiliar para estimacion de efecto de la flotación y su incertidumbre

#Factor de correccion por efecto de la flotación, expresado en término de densidades
Bu<-function(den_as,den_ac,den_c,den_s){
  (1-den_ac/den_c)/(1-den_as/den_s)
}

#Coeficientes de sensibilidad de Bu
dBu<-function(den_as,den_ac,den_c,den_s){
  .e1 <- 1 - den_as/den_s
  .e2 <- .e1^2
  .e3 <- 1 - den_ac/den_c
  c(den_as = .e3/(den_s * .e2), den_ac = -(1/(den_c * .e1)),
    den_c = den_ac/(den_c^2 * .e1), den_s = -(den_as * .e3/(den_s^2 * .e2)))
}

#===============================================================================
#  Curvas de titulación: Determinación de puntos finales
#===============================================================================


#' determina el punto final en una curva de titulación
#'
#' Obtiene el punto final de una curva de titulación en el intervalo expecificado, empleando diferentes enfoques:
#' i)   Máximo de la primera derivada numérica. La incertidumbre se obtiene asumiendo un distribución triangular entre los dos puntos adyacentes al valor encontrado, dentro del conjunto de datos crudos.
#' ii)  Punto de inflexión obtenido empleando estimador de distancia extrema del paquete inflection (1)
#' iii) Punto de inflexión empleando modelo de regresión logística empleando paquete nplr
#'
#' La incertidumbre devuelta corresponde al resultado de la estimación empleando los algoritmos indicados.
#'
#' @param x_ conjunto de datos x asociados a la curva de titulación
#' @param y_ conjunto de datos y asociados a la curva de titulación
#' @param intervalo_x corresponde al intervalo en el que busca el punto de inflexión de la curva dada

#'
#' @return data.frame con factor con abcisa en el punto de equivalencia e incertidumbre estándar asociada a la determinación
#'
#' @importFrom inflection edeci
#' @import nplr
#'
#' @seealso [edeci()],[nplr()]
#'
#' @references
#' (1) Demetris T. Christopoulos, On the Efficient Identification of an Inflection Point, International Journal of Mathematics and Scientific Computing, 6 (1), 13-20, 2016.
#' \url{https://veltech.edu.in/wp-content/uploads/2016/04/Paper-04-2016.pdf}.
#'
#' @export
puntoFinal<-function(x_,y_,intervalo_x=c(0,0)){

  # intervalo_x=c(5,15)
  # x_=datosAnalisis$valoracionNaOH$V
  # y_=datosAnalisis$valoracionNaOH$pH
  x=x_
  y=y_
  if(sum(intervalo_x)>0.0){
    .indx_inf=min(which(x>=intervalo_x[1] & x<=intervalo_x[2]))
    .indx_sup=max(which(x>=intervalo_x[1] & x<=intervalo_x[2]))
    x=x_[.indx_inf:.indx_sup]
    y=y_[.indx_inf:.indx_sup]
  }

  #========================================
  # PE empleando primera derivada numérica
  #=======
  df=forwardDerivative(x,y)
  y_range=max(y)-min(x)
  ind_max=which.max(df$dy_dx)
  PE_dnum=df$x_avg[ind_max]
  u_PE_dnum=(x[ind_max+1]-x[ind_max])/sqrt(6)
  y_pos_dnum=(y[ind_max]+y[ind_max+1])/2

  #===================================
  # PE empleando regresion logistica
  #======
  y.prop=convertToProp(y)
  model.logR=nplr(x,y.prop,useLog = FALSE,npars = 5)
  #plot(model.logR, ylim = range(0, 1))
  rNPLR=summary(model.logR)
  PE_nplr=as.numeric(rNPLR$value[13])
  y_nplr=as.numeric(rNPLR$value[14])
  #y_pos_nplr=y_nplr*y_range
  y_pos_nplr=(y[max(which(x<=PE_nplr))]+y[min(which(x>=PE_nplr))])/2
  u_68=getEstimates(model.logR,targets = c(y_nplr),conf.level = (  1-2*pnorm(-1,lower.tail = TRUE)  ))
  u_PE_nplr=(u_68$x.841-u_68$x.158)/2

  #================================================
  # PE empleando estimador de distancia extrema
  #======

  model.ede=edeci(x,y,index = 0,k = 1)
  PE_ede=model.ede[3]
  u_PE_ede=(model.ede[6]-model.ede[5])/2
  y_pos_ede=(y[max(which(x<=PE_ede))]+y[min(which(x>=PE_ede))])/2

 #Resultados
 plot(x,y,type="l",col="blue",main="Determinacion de punto de inflexion en curva de titulacion")
 points(x,y,pch=16,cex=0.5)
 fT=(max(y)-min(y))/(max(df$dy_dx)-min(df$dy_dx))
 lines(df$x_avg,df$dy_dx*fT+min(y),type="l",lty=1,col="gray")
 abline(v=c(PE_dnum,PE_nplr,PE_ede),lty=2,col=c("red","cyan","green"))


 points(PE_dnum,y_pos_dnum,pch=20,cex=.75,col="red")
 arrows(PE_dnum-u_PE_dnum,y_pos_dnum,PE_dnum+u_PE_dnum,y_pos_dnum, length=0.05,angle = 90, code=3, col="red")

 points(PE_nplr,y_pos_nplr,pch=20,cex=.75,col="cyan")
 arrows(PE_nplr-u_PE_nplr,y_pos_nplr,PE_nplr+u_PE_nplr,y_pos_nplr, length=0.05,angle = 90, code=3, col="cyan")

 points(PE_ede,y_pos_nplr,pch=20,cex=.75,col="green")
 arrows(PE_ede-u_PE_ede,y_pos_ede,PE_ede+u_PE_ede,y_pos_ede, length=0.05,angle = 90, code=3, col="green")

 legend("topleft",c(paste("Derivada numerica:",format(round(PE_dnum,2),nsmall = 2)),
                    paste("Regresion logistica: ",format(round(PE_nplr,2),nsmall = 2)),
                    paste("Estimador de maxima distancia: ",format(round(PE_ede,2),nsmall = 2))),
                    lty = 2,col=c("red","cyan","green"),bty ="n",cex=0.85)

  list("d1num"=data.frame(PE=PE_dnum,u_PE=u_PE_dnum),"rLog"=data.frame(PE=PE_nplr,u_PE=u_PE_nplr),
      "ede"=data.frame(PE=PE_ede,u_PE=u_PE_ede))

}

# Funcion auxiliar para funcion de punto final
forwardDerivative<-function(x,y){
  len=length(x)
  Dy=y[2:len]-y[1:len-1]
  Dx=x[2:len]-x[1:len-1]
  x_avg=(x[2:len]+x[1:len-1])/2
  dy_dx=Dy/Dx
  data.frame(dy_dx,x_avg,Dy,Dx)
}

#===============================================================================
#  Pureza de reactivos
#===============================================================================

#' Estimado de pureza e incertidumbre basados en la ficha técnica del estándar.
#'
#' Devuelve estimados de puereza en incertidumbre basados en especificaciones técnicas del reactivo (1-2)
#' Los valores se deben ingresar como fracciones másicas en el intervalo 0-1. Sí se trata de un MRC no es necesario emplear esta función.
#'
#' @param P pureza, pureza minima o vector con pureza minima y pureza máxima c(P_min,P_max), expresadas como fracción másica.
#' @param delta variación indicada sin especificación de factor de cobertura o tipo de distribución. Si delta=0, se asume que P es una pureza mínima
#' @param tecnica técnica de análisis descrita: "titulacion", "cromatografia", "n" (no indicada)
#' @param sInstrum indica si el delta proveido corresponde a la incertidumbre instrumental (FALSE, por omisión)
#' @param dist distribucion asumida para casos donde se especifica un estimado de P y delta sin conocimiento de tipo de distribucion o factor de cobertura; o valores máximos y mínimos de P. Por omisión se asume una distribución "rectangular".
#' Sí solo se indica un valor informativo estimado de P, sin información sobre incertidumbre, use delta=1E-10 y dist="triangular" para obtener un estimado de la incertidumbre estándar con valor centrado en P y distribución triangular en el intervalo +/-(1-P).
#'
#' @return data.frame con pureza e incertidumbre estimadas a partir de la especificaciones
#'
#' @examples
#' pureza(P=0.999) # Solo se indica el límite inferior de la pureza
#'
#' @references
#' (1) van Look, G., Meyer, V.R., The purity of laboratory chemicals with regard to measurement uncertainty, 	Analyst, 127, 825-829, 2002. \url{https://doi.org/10.1039/B107958C}
#'
#'  (2) Borges, R., Meyer, V.R. The uncertainty of purity of reference materials must be known, J. of Pharm. and Biomed. Anal., 77, 40-43, 2013. \url{https://doi.org/10.1016/j.jpba.2013.01.008}
#' @export
pureza<-function(P,delta=0.0,tecnica="n",sInstrum=FALSE,dist="rectangular"){
  len=length(P)
  .delta=abs(delta)
  if(len==1){
    if(.delta<=1E-15){
      if(tecnica=="titulacion"){
        .P=(1.01+P)/2
        .u_P=(1.01-P)/(2*sqrt(3))
      }else{
        .P=P+(1.0-P)/sqrt(2)
        .u_P=(1.0-P)/sqrt(18)
      }
    }else{
      if(sInstrum){
        .P=P
        .u_P=sqrt(sum(((1-P)/factoresEstadarizacion(dist))^2+delta^2))
      }else{
        .P=P
        .u_P=delta/factoresEstadarizacion(dist)
      }
    }
  }else{
    .P=(P[1]+P[2])/2
    .u_P=(P[1]-P[2])/(2*factoresEstadarizacion(dist))
  }
}




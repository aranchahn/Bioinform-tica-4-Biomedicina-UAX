### TRABAJO ORDINARIA BIOINFORMÁTICA R ###
> # 1.
> datos <- read.delim("/Users/ariihervaas/Desktop/Trabajo FINAL de R-20261023/archivo-entrada.txt", header = TRUE, check.names = FALSE)
> head(datos)
  Gene Name    C2-TPM    C1-TPM    C3-TPM   1W5-TPM   1W6-TPM   1W5-COV
1      urtA  30.40059  34.71058  32.35147  0.000000  0.000000  0.000000
2      recQ 101.01847 107.75805 107.83519  3.457582 10.163883  3.312236
3      kdsA  60.78179  57.62743  59.47843  0.000000 12.250913  0.000000
4      lysA 129.41655 109.94573 126.75361 24.608057 32.755714 23.573610
5       sir 153.05125 154.52052 160.49197  2.828362  8.774453  2.709466
6       map  51.67744  54.92715  55.71499  0.000000  5.603065  0.000000
> tail(datos)
    Gene Name     C2-TPM   C1-TPM    C3-TPM    1W5-TPM    1W6-TPM    1W5-COV
614      holB  26.528334 25.11444  25.84524   0.000000  14.037302   0.000000
615      xseA   4.958211 10.49672  13.32347   0.000000   8.579851   0.000000
616      ccsB  79.833740 52.41565  56.37512 108.076942  77.530563 103.533722
617      glnA 100.763817 97.39648  90.88831  78.461060 109.376648  75.162804
618      alaS  27.634129 23.34384  26.04073   9.950721  15.988743   9.532423
619      hisH 145.615585 68.18367 105.57125  44.589828 103.541252  42.715408
> dim(datos)
[1] 619   7

> # 2.
# Media y Desviación Estándar para elControl (columnas 2 a 4)
> datos$media_control <- rowMeans(datos[, 2:4])
> datos$DE_control <- apply(datos[, 2:4], 1, sd)
# Media y Desviación Estándar para la Condición (columnas 5 a 7)
> datos$media_condicion <- rowMeans(datos[, 5:7])
> datos$DE_condicion <- apply(datos[, 5:7], 1, sd)
# Comprobar como queda la tabla con las 4 columnas
> head(datos)
  Gene Name    C2-TPM    C1-TPM    C3-TPM   1W5-TPM   1W6-TPM   1W5-COV media_control DE_control
1      urtA  30.40059  34.71058  32.35147  0.000000  0.000000  0.000000      32.48754   2.158215
2      recQ 101.01847 107.75805 107.83519  3.457582 10.163883  3.312236     105.53724   3.913556
3      kdsA  60.78179  57.62743  59.47843  0.000000 12.250913  0.000000      59.29588   1.585084
4      lysA 129.41655 109.94573 126.75361 24.608057 32.755714 23.573610     122.03863  10.557057
5       sir 153.05125 154.52052 160.49197  2.828362  8.774453  2.709466     156.02125   3.940839
6       map  51.67744  54.92715  55.71499  0.000000  5.603065  0.000000      54.10653   2.140215
  media_condicion DE_condicion
1        0.000000     0.000000
2        5.644567     3.914517
3        4.083638     7.073068
4       26.979127     5.029338
5        4.770760     3.467809
6        1.867688     3.234931
> # 3.
> png("03_histogramas.png", width = 1200, height = 600)
> 
> par(mfrow = c(1, 2))
> hist(datos$media_control, col = "skyblue", main = "Control",
+      xlab = "Media de expresion", ylab = "Numero de genes")
> hist(datos$media_condicion, col = "salmon", main = "Condicion",
+      xlab = "Media de expresion", ylab = "Numero de genes")
> dev.off()
null device 
          1 
> par(mfrow = c(1, 2))
> 
> hist(datos$media_control, col = "skyblue",
+      main = "Control",
+      xlab = "Media de expresión", ylab = "Número de genes")
> 
> hist(datos$media_condicion, col = "salmon",
+      main = "Condición",
+      xlab = "Media de expresión", ylab = "Número de genes")

> # 4. 
> png("04_boxplot.png", width = 1000, height = 700)
> boxplot(datos$media_control, datos$media_condicion,
+ names = c("Control", "Condicion"), col = c("skyblue", "salmon"),
+ main = "Expresion media de los genes", ylab = "Media de expresion")
> legend("topleft", legend = c("Control", "Condicion"),
+        fill = c("skyblue", "salmon"), bty = "n")
> dev.off()
null device 
          1 
> par(mfrow = c(1, 1))
> 
> boxplot(datos$media_control, datos$media_condicion,
+         names = c("Control", "Condición"),
+         col = c("skyblue", "salmon"),
+         main = "Expresión media de los genes",
+         ylab = "Media de expresión")
> 
> legend("topleft",
+        legend = c("Control", "Condición"),
+        fill = c("skyblue", "salmon"),
+        bty = "n")

> # 5.
> shapiro.test(datos$media_control)

	Shapiro-Wilk normality test

data:  datos$media_control
W = 0.06848, p-value < 2.2e-16

> shapiro.test(datos$media_condicion)

	Shapiro-Wilk normality test

data:  datos$media_condicion
W = 0.032856, p-value < 2.2e-16

> datos$p_value <- NA
> for (i in 1:nrow(datos)) {
+   control <- as.numeric(datos[i, 2:4])
+   condicion <- as.numeric(datos[i, 5:7])
+   # P exacto si no hay valores repetidos; aproximado si hay empates.
+   sin_empates <- length(unique(c(control, condicion))) == 6
+   datos$p_value[i] <- wilcox.test(control, condicion,
+                                 exact = sin_empates)$p.value
+ }
> head(datos)
  Gene Name    C2-TPM    C1-TPM    C3-TPM   1W5-TPM   1W6-TPM   1W5-COV media_control DE_control
1      urtA  30.40059  34.71058  32.35147  0.000000  0.000000  0.000000      32.48754   2.158215
2      recQ 101.01847 107.75805 107.83519  3.457582 10.163883  3.312236     105.53724   3.913556
3      kdsA  60.78179  57.62743  59.47843  0.000000 12.250913  0.000000      59.29588   1.585084
4      lysA 129.41655 109.94573 126.75361 24.608057 32.755714 23.573610     122.03863  10.557057
5       sir 153.05125 154.52052 160.49197  2.828362  8.774453  2.709466     156.02125   3.940839
6       map  51.67744  54.92715  55.71499  0.000000  5.603065  0.000000      54.10653   2.140215
  media_condicion DE_condicion    p_value
1        0.000000     0.000000 0.06360257
2        5.644567     3.914517 0.10000000
3        4.083638     7.073068 0.07652250
4       26.979127     5.029338 0.10000000
5        4.770760     3.467809 0.10000000
6        1.867688     3.234931 0.07652250

> # 6. 
> genes_DE <- datos[datos$p_value < 0.05, ]
> genes_DE[, c("Gene Name", "media_control", "media_condicion", "p_value")]
[1] Gene Name       media_control   media_condicion p_value        
<0 rows> (o 0- extensión row.names)
> nrow(genes_DE)
[1] 0
> sobreexpresados <- genes_DE[genes_DE$media_condicion > genes_DE$media_control, ]
> reprimidos <- genes_DE[genes_DE$media_condicion < genes_DE$media_control, ]
> nrow(sobreexpresados)
[1] 0
> nrow(reprimidos)
[1] 0
> # No se han encontrado genes con p_value < 0.005. Solo hay tres valores por grupo y el test Wilcoxon bilateral. Tiene poca capacidad para detectar diferencias con tan pocos datos, pero no sognifica que existan cambios de expresión, sino que se detectan como significativos con este análisis.

> # 7.
>filtered <- as.matrix(genes_DE[, 2:7]) 
rownames(filtered) <- genes_DE$`Gene Name`
rowv <- as.dendrogram(hclust(as.dist(1 - cor(t(filtered))))) 
colv <- as.dendrogram(hclust(as.dist(1 - cor(filtered)))) 
heatmap(filtered, Rowv = rowv, Colv = colv, cexCol = 0.7, labRow = FALSE) 
if (!require("gplots")) install.packages("gplots") 
if(!require("RColorBrewer")) install.packages("RColorBrewer") 
library(gplots) 
library(RColorBrewer) 
heatmap.2(filtered, Rowv = rowv, Colv = colv, cexCol = 0.7, col = rev(redblue(256)), scale = "row", labRow = FALSE) 
>pdf("GSE5583_DE_Heatmap.pdf") 
heatmap.2(filtered, Rowv = rowv, Colv = colv, cexCol = 0.7, col = rev(redblue(256)), scale = "row", labRow = FALSE) 
dev.off() 
>write.table(filtered, "GSE5583_DE.txt", sep = "\t", quote = FALSE)
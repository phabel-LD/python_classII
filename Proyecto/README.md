---
title: "README.md"
author: "Phabel Antonio López Delgado & Daianna González Padilla"
date: "12/1/2021"
---

# Proyecto BioPythonII

## Descripción
En este repositorio se encuentra todos los archivos y scripts del proyecto semestral, así como el reporte del mismo.

## Requerimientos
Python v3
R v3.4.0

## Licencia
Para ver la licencia ver al archivo [LICENCIA](LICENCIA)

## Contacto
Phabel Antonio López Delgado <phabel@lcg.unam.mx>
Daianna González Padilla <daianna@lcg.unam.mx>

## Reporte

### Métodos
Los métodos utilizados se enfocan en un análisis funcional y estructural de un conjunto de proteínas a través de una comparación de sus motivos.

El primer paso es la descarga de archivos .pdb con la información de las proteínas a ser analizadas. Posteriormente se realiza una búsqueda de motivos con relevancia biológica para cada proteína analizada.

Una vez con los datos de tipos de motivos y sus respectivas frecuencias contenidas en las proteínas de interés, se pide al usuario ingresar el nombre de una proteína (str) como referente para el análisis. Para ello, se pueden utilizar argumentos adicionales para el análisis con el fin de utilizar librerías adicionales para mostrar los resultados (numpy, seaborn, matplotlib.pyplot, pandas).

1) Obtención y cuantificación de motivos
EL primer punto es cuantificar los motivos de interés biológico, particularmente puentes disúlfido (S-S), hélices alfa y hojas beta plegadas, a partir de la información particular que se puede hallar en los archivos .pdb. También se pueden buscar motivos personalizados a través de un parseo con expresiones regulares (regex). La forma de guardar toda la información encontrada para su posterior análisis es en una lista de diccionarios, donde cada diccionario almacena la información extraída de una proteína; dando un diccionario por proteína. La llaves del diccionario representarán los nombres de los motivos y los valores sus frecuencias halladas.

2) Consecutivamente, se pedi el nombre (ID) (str) de una proteína query a ser analizada funcional y estructuralmente a partir de su frecuencia y composición de motivos. Dada una proteína query, su diccionario con el análisis de sus motivos será comparado con el resto de los diccionarios que representan a las otras proteínas. De esta forma se puede determinar semejanza y relaciones estructurales, funcionales y por lo tanto, evolutivas.

3) El siguiente punto es generalizar el análisis de intersección pairwise, donde el modelo de análisis propuesto compara motivo a motivo, la intersección relativa de las proteínas. Esto se hace basándose en la lógica del índice de Jaccard para comparar conjuntos, pero con un ligero ajuste: se dividirá el valor mínimo entre el máximo para cada motivo. Esta generalización comparara la proteína query contra todas las proteínas obtenidas de los .pdb. Dando como resultado una matrix de intersecciones.

4) Para conseguir más información y un score útil para poder comparar directamente la semejanze estructural, funcional y evolutiva, se obtendrá un vector de promedios iterando sobre las filas de la matriz de intersección.

5) Análisis opcionales que pueden reaizarse proveyéndo los argumentos necesarios, incluyen presentar la matriz de intersección y el vector de promedios como dataframes, plotear un heatmap a partir de la matriz de intersección, y obtener la proteína con el mejor score (promedio) a partir del vector de promedios y presentarla como el best match para la proteína query de entre toda la lista de proteínas.

6) Nótese que con el fin de evitar la autocomparación, el score de la proteína query consigo misma se trivialia a cero.

### Pregunta
Este programa pretende proporcionar un herramienta rápida y simple para la identificación y conteo de motivos de una secuencia proteíca, y la subsecuente compración 
de una proteína query con otras con base en el número de cada uno de los motivos econtrados. La pregunta que se pretende responder es si existe una relación o 
similitud funcional entre las proteínas que tienen numeros parecidos o cercanos de cada uno de los motivos que poseen. Para ello se deberán tomar casos conocidos de
proteínas funcionalmente parecidas y usando este programa calcular su número de motivos para verificar si poseen o no numeros parecidos.   
 

### Resultados
Dependiendo de los argumentos de análisis pasados por línea de comandos al programa, el número de resultados obtenidos puede variar.

Con la línea de ejemplo base...

python3 ProteinAnalysis.py -i 1kcw.pdb 1fat.pdb 3jbz.pdb 1kbe.pdb 4g68.pdb 5v6r.pdb -d 2 -a default -b default --int_matx_df --means_df --best --heatmap

Se obtienen los motivos...
[Análisis de motivos](Resultados/motifs.png)


Ingresando la proteína query...
[Proteína query ingresada](Resultados/query_prot.png)


Se crea la matriz de intersección como objeto numpy por default...
[Matriz de intersección como array de numpy](Resultados/matrix_np.png)


Se crea el vector de promedios como objeto numpy por default...
[Vector de promedios como array de numpy](Resultados/vector_np.png)


Se puede convertir la matriz de intersección de array de numpy a dataframe de pandas (opcional --int_matx_df)...
[Matriz de intersección como dataframe de pandas](Resultados/matrix_df.png)


Se puede convertir el vector de promedios de array de numpy a dataframe de pandas (opcional --means_df)...
[Vector de promedios como dataframe de pandas](Resultados/vector_df.png)


Se puede mostrar en pantalla el best match para la proteína query tomando el máximo score obtenido en el vector de promedios (opcional --best)...
[Best match encontrado](Resultados/best.png)


Se puede imprimir el heatmap correspondiente a la matriz de intersección (opcional --heatmap)...
[Heatmap como objeto de seaborn](Resultados/heatmap.png)


Dando por finalizado los resultados que se pueden obtener. Los análisis opcionales se pasan con argumentos adicionales; véase en la documentación del script.

### Conclusión
De la ejecución del programa es posible concluir que el próposito del proyecto se cumplió al crear una herramienta rápida, simple y eficiente para identificar y contar
motivos, así como para comprar proteínas según tal conteo. Sin embargo, la hipótesis de que existe una relación funcional según la similitud en número de motivos 
no es clara ni contundente puesto que el reconocimiento mismo de motivos se complica al no existir patrones únicos o consenso para la identificación de los mismos; se necesitarán estudios evolutivos y funcionales que integren además información sobre la localización de los motivos y la similitud en sus secuencias de aminoácidos.



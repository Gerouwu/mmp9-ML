# MMP-9: clasificación de bioactividad

Proyecto del taller de bioinformática: **SMILES TF-IDF y descriptores PaDEL**, con
**XGBoost, Random Forest y regresión logística**, 100 repeticiones por combinación.
Positivo = molécula activa frente a MMP-9 humana (CHEMBL321 / P14780).

## Consultar los resultados

- [Informe sencillo con figuras](reports/informe_rendimiento.html) (abrir en navegador).
- [Informe en Markdown](reports/informe_rendimiento.md), con referencias a los `.py`.
- [Tabla resumen](reports/rendimiento_resumen.csv).
- [Todas las imágenes PNG en ZIP](reports/imagenes_mmp9.zip): matrices de las seis
  combinaciones en conteos medios y porcentajes, gráficas individuales y comparaciones.
  Regenerar con `python -m mmp9.plots`; detalles en `reports/imagenes/README.md`.
- `results/main/`: métricas por ejecución, predicciones, ajuste de hiperparámetros,
  particiones, versiones, modelos de la semilla 0 y auditoría de resultados.

El informe IEEE queda para una fase posterior. El informe presente describe los
resultados computacionales y sus limitaciones; no demuestra eficacia farmacológica.

## Organización

```text
mmp9/
  data.py         Descarga ChEMBL y curación trazable de las mediciones
  features.py     Descriptores PaDEL 1D/2D (requiere Java)
  experiment.py   Selección interna y 600 evaluaciones externas
  report.py       Informe HTML/Markdown, tablas y figuras
  verify.py       Auditoría independiente de los resultados guardados
data/
  raw/            Copia de respuestas de ChEMBL, objetivo y procedencia
  processed/      Moléculas, exclusiones, conflictos y descriptores
results/main/     Resultados del experimento completo
reports/          Informe sencillo y figuras
tests/            Pruebas de comportamiento
```

La imagen del taller y el `.py` docente permanecen intactos en la raíz. Ese `.py`
es una exportación de Colab, se conserva como referencia y no se ejecuta.
Los entornos, Java descargado y cachés no se incluyen en Git.

## Reproducir

Requisitos: Python 3.12 y Java 8+ en PATH. Para esta ejecución se usó Temurin JRE 8,
localizado en `.tools/java`; `features.py` también reconoce esa ubicación en Windows.
Java puede obtenerse de [Eclipse Adoptium](https://adoptium.net/temurin/releases/).

Desde la raíz, en PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pytest -q
# Usa los datos originales incluidos (no requiere nueva descarga):
.\.venv\Scripts\python.exe -m mmp9.data
# Opcional si se quiere recalcular PaDEL; los descriptores ya están incluidos:
.\.venv\Scripts\python.exe -m mmp9.features
.\.venv\Scripts\python.exe -m mmp9.experiment --repeats 100 --workers 3 --output results/reproduccion
.\.venv\Scripts\python.exe -m mmp9.verify --output results/reproduccion
.\.venv\Scripts\python.exe -m mmp9.report --results results/reproduccion
```

Alternativa integrada, que aprovecha los descriptores incluidos si corresponden
a los datos: `python run_pipeline.py --output results/reproduccion`. Esta instrucción
ejecuta curación, características (cuando hacen falta), entrenamiento, auditoría e informe.

Para actualizar desde ChEMBL: `python -m mmp9.data --download`. Esto reemplaza la
copia de datos y exige recalcular descriptores. Una nueva versión de ChEMBL puede
cambiar los resultados. Para reproducir los resultados entregados usar la copia
incluida. La generación del informe reemplaza los archivos de `reports/`.

Prueba corta: `python -m mmp9.experiment --repeats 2 --output results/smoke`.
Cada semilla completada crea un punto de reanudación local (ignorado por Git).
Se rechaza reanudar si cambian datos, código, versiones o configuración; usar otra
carpeta de salida. Los resultados consolidados y particiones sí están versionados.

## Decisiones del experimento

- IC50 exacto en nM, positivo y finito; ensayos de unión con confianza 9,
  sin alertas de validez ni duplicados potenciales.
- Estructura padre canónica, mediana de mediciones por estructura y exclusión de
  moléculas con mediciones tanto ≤1000 como ≥10000 nM.
- Activa: IC50 ≤1000 nM. Inactiva: IC50 ≥10000 nM. Las intermedias se conservan
  en el conjunto curado, pero se excluyen del experimento binario.
- División estratificada externa 80/20, semillas 0–99. División interna 75/25
  del entrenamiento para elegir entre dos configuraciones por modelo usando F1.
  Son 1200 ajustes internos y 600 reajustes/evaluaciones externas.
- Vocabulario, imputación, selección por varianza y escalado se ajustan solo sobre
  entrenamiento; se reajustan con el 80% tras seleccionar hiperparámetros.
- Todas las combinaciones comparten moléculas y particiones. Umbral fijo 0.5.
- Media y desviación estándar muestral de accuracy, F1, sensibilidad y especificidad.
- No hay separación por familia química ni validación externa; el informe discute
  el posible optimismo de la partición aleatoria y el solapamiento entre repeticiones.

Los modelos `.joblib` incluyen preprocesador y clasificador, corresponden a semilla
0 y deben cargarse solo desde una fuente confiable. XGBoost incluye además `.ubj`;
el archivo nativo por sí solo no incluye la transformación de entradas.

## Fuentes y atribución

Datos de [EMBL-EBI ChEMBL](https://www.ebi.ac.uk/chembl/), licencia CC BY-SA 3.0.
Objetivo y versión verificados en `data/raw/target.json` y `chembl_status.json`;
fecha, consultas y hashes en `provenance.json`. Las tablas derivadas mantienen la
atribución y licencia de los datos de origen. Código docente original atribuido
a los autores indicados en su encabezado; se conserva sin modificaciones.
Documentación de las bibliotecas enlazada en el informe.

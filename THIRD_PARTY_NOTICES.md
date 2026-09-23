# Licencias y atribuciones de terceros

La licencia MIT de este repositorio se aplica al código y a la documentación
originales del proyecto. No sustituye las licencias ni los derechos de los
materiales de terceros descritos a continuación.

## Datos de ChEMBL

Los datos proceden de EMBL-EBI ChEMBL, versión 37, objetivo CHEMBL321 (MMP-9
humana). Los datos originales y sus tablas derivadas se distribuyen bajo
[Creative Commons Attribution-ShareAlike 3.0 Unported](https://creativecommons.org/licenses/by-sa/3.0/),
no bajo MIT. Véase la [información oficial de licencia de ChEMBL](https://chembl.gitbook.io/chembl-interface-documentation/about).

Los archivos de `data/raw/` conservan los registros y la procedencia de la
consulta. La preparación de `data/processed/` filtra mediciones IC50, normaliza
estructuras, agrupa mediciones, calcula descriptores y asigna clases. Las
predicciones y tablas de resultados que incorporan estos datos deben conservar
su atribución y las condiciones aplicables de ChEMBL. La metodología y las
transformaciones se documentan en el README y en los informes.

## Material docente de referencia

`mib_introduction_drug_discovery_through_ml.py` es material docente de referencia
conservado sin modificaciones. Su encabezado atribuye la autoría a Alexander
Rodríguez López y Alvaro David Orjuela Cañón, e indica una adaptación del
[material de Data Professor](https://github.com/dataprofessor/bioinformatics_freecodecamp/blob/main/CDD_ML_Part_1_Acetylcholinesterase_Bioactivity_Data_Concised.ipynb).
La licencia MIT del proyecto no se otorga sobre este archivo ni sobre la imagen
original del enunciado `WhatsApp Image 2026-09-22 at 6.00.26 PM.jpeg`.
Sus autores y titulares conservan sus derechos; su inclusión no implica una
autorización adicional de reutilización ni una cesión de marcas o logotipos.

## Dependencias

Las bibliotecas y herramientas utilizadas conservan sus respectivas licencias.
La licencia MIT del proyecto no modifica las condiciones de XGBoost,
scikit-learn, RDKit, PaDEL-Descriptor ni de las demás dependencias.

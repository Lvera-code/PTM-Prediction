# PTM-Prediction

[![tests](https://github.com/Lvera-code/PTM-Prediction/actions/workflows/tests.yml/badge.svg)](https://github.com/Lvera-code/PTM-Prediction/actions/workflows/tests.yml)

Pipeline de prediccion de zonas de modificacion post-traduccional (PTM) a
partir de FASTA o PDB/mmCIF. Segundo proyecto del CNB (encargado por Carlos
Oscar Sorzano), independiente de
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction).

## Arquitectura

Numeracion de fases alineada con el proyecto 1
([BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)):
Fase 1 -> Fase 1.5 (PDB unicamente) -> Fase 2 (motores) -> Fase 3 (nucleo) -> Fase 3b
(cruces informativos, PDB unicamente).

- **Camino FASTA**: Fase 1 (saneamiento) -> Fase 2: **DeepMVP** (motor unico,
  6 tipos de PTM: fosforilacion, acetilacion, metilacion, sumoilacion,
  ubiquitinacion, N-glicosilacion) -> Fase 3.
- **Camino PDB**: Fase 1.5 (extraccion de secuencia ATMSEQ + mapeo de
  posiciones via `gemmi`) -> Fase 2: consenso **DeepMVP + DeepPTMPred** (17
  tipos de PTM; DeepPTMPred exige `pdb_path` obligatorio, sin modo
  solo-secuencia, de ahi la asimetria entre caminos) -> Fase 3 -> **Fase 3b**.
- **Fase 3 (nucleo)**: anotacion/filtrado + logica de decision de flujo sobre
  las predicciones crudas de Fase 2 (`src/engines/ptm_annotation.py`).
  Fusiona consenso donde DeepMVP y DeepPTMPred coinciden en tipo+posicion.
  `n_linked_glycosylation` esta deliberadamente EXCLUIDO de ESA fusion
  especifica (modelo de DeepPTMPred sin poder discriminativo real para este
  tipo, AUROC 0.51 -- ver STATUS.md), aunque ambos motores lo siguen
  reportando por separado -- en su lugar, Camino PDB tiene un consenso REAL
  distinto para este tipo especifico: **DeepMVP + EMNGly + StackGlyEmbed**
  (decision 2026-08-06, reemplaza al CoNglyPred original, sin pesos
  publicados -- ver STATUS.md "Decision 2"), `pasa_umbral` = al menos 1 de
  los motores disponibles pasa, `consenso` = al menos 2 de 3.
- **Corroboracion opcional** (nunca deciden `pasa_umbral`/`consenso`, solo
  informan sobre sitios que el consenso YA acepto): **MeToken** (tipo, Camino
  PDB), **StackGlyEmbed** (N-glicosilacion, Camino FASTA -- en Camino PDB es
  motor de consenso, ver arriba), **evidencia de via secretora via UniProt**
  (columna `via_secretora_evidencia`, N-glicosilacion, ambos caminos) y
  **aviso de competencia/crosstalk entre PTMs** (columna
  `ptm_crosstalk_aviso`, cuando 2+ tipos que compiten por el mismo grupo
  quimico de un residuo pasan en la misma posicion) -- ver "Alcance e
  interpretacion" abajo.

## Alcance e interpretacion

Analisis de coherencia biologica 2026-08-07 (equivalente, para este
proyecto, a las recomendaciones de Carmen Elena Gómez para
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)):
huecos de comunicacion de alcance encontrados, ninguno un bug de codigo.

- **Capacidad predicha, no ocurrencia observada.** Ningun motor de este
  pipeline (DeepMVP/DeepPTMPred/EMNGly/StackGlyEmbed) modela la via
  biosintetica real del sustrato -- todos predicen desde secuencia/
  estructura si un sitio ES modificable en principio, nunca si esa PTM
  ocurre realmente en una celula/tejido/condicion especifica (eso depende
  de que la enzima real este co-expresada y co-localizada con el sustrato,
  algo que ningun predictor de secuencia puede capturar). Un resultado de
  este pipeline debe leerse como "sitio potencialmente modificable", no
  como "esta proteina SE modifica". Impreso como aviso al final de cada
  corrida (`pipeline.py::INTERPRETATION_DISCLAIMER`).
- **Sin filtro de via secretora para N-glicosilacion.** La N-glicosilacion
  ocurre quimicamente en el lumen del RE/Golgi -- una proteina puramente
  citoplasmatica/nuclear practicamente nunca se glicosila en un sequon,
  sin importar el score de consenso. La columna `via_secretora_evidencia`
  (via `src/structural/uniprot_localization_client.py`, evidencia REAL de
  UniProt, nunca inventada) avisa cuando el consenso acepto un sitio
  N-glico sin evidencia conocida de localizacion secretora -- NUNCA
  rechaza el candidato (mismo criterio que la recomendacion de Carmen
  Elena de no descartar candidatos de N-glico en el proyecto 1). Alcance
  limitado a `n_linked_glycosylation` -- `o_linked_glycosylation` tiene dos
  vias biologicas distintas (O-GlcNAc citoplasmatica/nuclear vs mucina en
  la via secretora) que este cliente no distingue.
- **Sin modelado de competencia entre PTMs del mismo residuo.** Varios
  tipos modifican el mismo grupo quimico de un residuo y son mutuamente
  excluyentes en una misma molecula/instante (acilo-lisina: acetilacion/
  ubiquitinacion/sumoilacion/metilacion/malonilacion/glutarilacion/
  succinilacion/crotonilacion; tiol de cisteina: S-nitrosilacion/
  glutationilacion; guanidino de arginina: metilacion/citrulinacion;
  hidroxilo de Ser/Thr: fosforilacion/O-glicosilacion, la hipotesis
  "Yin-Yang"). La columna `ptm_crosstalk_aviso` avisa cuando 2+ tipos en
  competencia real pasan en la misma posicion -- ver
  `src/engines/ptm_annotation.py::_PTM_COMPETITION_GROUPS`.
- **Especificidad de quinasa para fosforilacion.** Ni DeepMVP ni DeepPTMPred
  distinguen QUE familia de quinasa fosforila un sitio -- ambos predicen
  "fosforilable en general". A diferencia del tipo de cadena de
  poliubiquitina (K48 = degradacion proteasomal; K63 = senalizacion/
  reparacion, NO degradativo -- depende de que E2/E3 real conjuga ubiquitinas
  adicionales dentro de la celula, un evento celular posterior no observable
  en secuencia/estructura), la especificidad de quinasa SI es una propiedad
  local de secuencia -- existe una fuente real y publicada:
  [Kinase Library](https://kinase-library.phosphosite.org) (Johnson et al.
  2023 *Nature*, 303 quinasas Ser/Thr + Yaron-Barir et al. 2024 *Nature*,
  kinoma Tyr completo). Para cada sitio de fosforilacion que el consenso
  acepta, las columnas `kinase_library_top_kinase`/
  `kinase_library_top_family`/`kinase_library_percentile`/
  `kinase_library_top3_kinases` (puramente informativas, nunca deciden
  `pasa_umbral`/`consenso`) reportan la quinasa/familia mas probable segun
  las matrices de especificidad publicadas. Verificado 2026-08-07 contra un
  sitio real (p53 S33): top hit ATM, coincide con la literatura real de
  respuesta a dano en el ADN.

## Estado actual (2026-08-10)

Pipeline completo end-to-end (Fase 1/1.5/2/3/3b) implementado, instalado y
verificado con corridas REALES (no solo planeado) en esta maquina: DeepMVP
(pesos + calibracion real), DeepPTMPred (PyRosetta + ESM-2 + calibracion real
de los 17 tipos, 2 bugs de train/inferencia encontrados y corregidos --
phi/psi y plDDT, ver STATUS.md), MeToken y StackGlyEmbed. Confirmado
feature-complete y sin Fase A/3c (modelado estructural real via PyRosetta,
eliminada del alcance) por feedback de Carlos -- ver STATUS.md.

**Decision 2 (segundo motor de consenso para N-glicosilacion) CERRADA**:
CoNglyPred (candidato original) confirmado sin pesos publicados en ningun
sitio -- reemplazado por **EMNGly** (pesos reales, verificados a nivel de
bytes), implementado y wireado al consenso de `n_linked_glycosylation` junto
con StackGlyEmbed (promovido de corroboracion informativa a motor de
consenso en Camino PDB). Ver STATUS.md "Decision 2" para el detalle
completo, incluyendo los 2 go/no-go checks (ambos PASARON, 2026-08-07): el
alineamiento de `structure_emb` se verifico contra sitios reales de GlyGen
en un PDB con huecos (Alpha-1-antitrypsin, 1QLP), y el MCC publicado
(0.736) se reprodujo -- de hecho se supero, 0.8197 -- corriendo el pipeline
real sobre el set independiente de N-GlyDE. El umbral
`EMNGLY_MIN_PROBABILITY=0.5` ya no es provisional.

## Uso

```bash
python pipeline.py --input inputs/mi_proteina.fasta
python pipeline.py --input inputs/mi_estructura.pdb   # Camino PDB: motores + consenso + cruces informativos
```

## Instalacion

```bash
pip install -r requirements.txt
pytest tests/
```

DeepMVP y DeepPTMPred requieren venvs/conda envs dedicados aparte (stacks
incompatibles entre si), nunca el venv principal del pipeline:

```bash
git clone https://github.com/bzhanglab/DeepMVP
conda create -n deepmvp python=3.7.10 -y
conda run -n deepmvp pip install -r DeepMVP/requirements.txt
export DEEPMVP_PYTHON_BIN=$(conda run -n deepmvp which python)
# Pesos: descarga manual desde https://deepmvp.ptmax.org/ -> DeepMVP/models/

git clone https://github.com/kuikui-wang/DeepPTMPred
conda create -n deepptmpred python=3.10 -y
# Instalar segun DeepPTMPred/pred/train_PTM/environment.yml (TensorFlow 2.15 +
# tensorflow-addons + torch + fair-esm) + pyrosetta-installer (ver README del
# repo) + checkpoint ESM-2 (~2.5GB, mas su companero *-contact-regression.pt)
export DEEPPTMPRED_PYTHON_BIN=$(conda run -n deepptmpred which python)
```

MeToken (corroboracion OPCIONAL e informativa del tipo, Camino PDB, ver
STATUS.md 2026-08-01 -- NUNCA un motor de consenso, se puede omitir sin
perder funcionalidad del pipeline principal, `Settings.METOKEN_ENABLED`
degrada solo con un aviso si no esta instalado):

```bash
git clone https://github.com/A4Bio/MeToken
conda create -n metoken python=3.10 -y
conda run -n metoken pip install --index-url https://download.pytorch.org/whl/cpu torch
conda run -n metoken pip install numpy scipy biopython transformers omegaconf tqdm pandas huggingface-hub h5py
conda run -n metoken pip install --no-build-isolation torch_scatter  # sin wheel prebuilt, compila desde fuente
export METOKEN_PYTHON_BIN=$(conda run -n metoken which python)
# Pesos: descarga release 1.0 (~88MB) y descomprime en MeToken/pretrained_model/
curl -L -o /tmp/pretrained_model.zip https://github.com/A4Bio/MeToken/releases/download/1.0/pretrained_model.zip
python -c "import zipfile; zipfile.ZipFile('/tmp/pretrained_model.zip').extractall('MeToken')"
```

StackGlyEmbed (N-glicosilacion, AMBOS caminos -- FASTA y PDB). En Camino
FASTA sigue siendo puramente informativa (`Settings.STACKGLYEMBED_ENABLED`
degrada solo con un aviso si no esta disponible). En Camino PDB fue
PROMOVIDA 2026-08-06 a motor real de consenso junto con EMNGly (ver abajo),
reemplazando el rol que hubiera tenido DeepPTMPred si no estuviera
confirmado muerto para `n_linked_glycosylation` -- sigue siendo opcional
(el pipeline no se cae sin ella), pero deja de ser meramente informativa
ahi. A diferencia de los demas motores, NO requiere instalacion propia
aqui -- reusa el venv/pickles YA instalados en el proyecto hermano
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)
(decision 2026-07-26: nunca se importa codigo entre proyectos, pero SI se
reusan recursos externos pesados ya instalados como venvs/pesos):

```bash
# Requiere que B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed
# ya exista (ver su propio README, Seccion 11). Si vive en otra ruta:
export STACKGLYEMBED_PYTHON_BIN=/ruta/a/B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed/bin/python
export STACKGLYEMBED_MODELS_DIR=/ruta/a/B-Cell-Epitope-Prediction/StackGlyEmbed/prediction
```

EMNGly (motor real de consenso para `n_linked_glycosylation`, Camino PDB
unicamente -- reemplaza a CoNglyPred, decision 2026-08-06, ver STATUS.md
"Decision 2": CoNglyPred confirmado sin pesos publicados en ningun sitio).
Opcional (`Settings.EMNGLY_ENABLED` degrada solo con un aviso, el consenso
de N-glicosilacion cae a DeepMVP+StackGlyEmbed si no esta disponible):

```bash
git clone https://github.com/StellaHxy/EMNgly
python3 -m venv .venv-emngly
.venv-emngly/bin/pip install fair-esm torch "scikit-learn==1.1.1" scipy pandas numpy tqdm wget
# 'wget' (paquete pip, no el binario CLI): dependencia transitiva real de
# MIF/sequence_models/trRosetta_utils.py (importado en cadena por
# pretrained.py aunque este proyecto nunca usa esa ruta de codigo) -- sin
# ella el import de MIF falla con ModuleNotFoundError (confirmado real
# 2026-08-07, no documentado en el environment.yml de EMNgly).
# pip resuelve numpy a la ultima 2.x por defecto, incompatible en runtime con
# el wheel compilado de scikit-learn==1.1.1 (ValueError: numpy.dtype size
# changed -- confirmado real 2026-08-07). Reinstalar despues, pinneado:
.venv-emngly/bin/pip install "numpy==1.23.5"
export EMNGLY_PYTHON_BIN=$(pwd)/.venv-emngly/bin/python

# MIF (embedding estructural, mif.pt) ya viene bundled en el clon -- sin
# descarga aparte.

# ESM-1b (~7.4GB) + companero de regresion de contactos (fair-esm SIEMPRE lo
# busca, mismo mecanismo que el checkpoint ESM-2 de DeepPTMPred arriba):
mkdir -p EMNgly/esm/checkpoints
curl -L -o EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S.pt \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
curl -L -o EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt \
  https://dl.fbaipublicfiles.com/fair-esm/regression/esm1b_t33_650M_UR50S-contact-regression.pt

# SVM (verificado 2026-08-06 con HTTP Range real, ver STATUS.md): elegido
# N-GlyDE.pickle (negativos restringidos al sequon, el regimen correcto),
# no N-GlyAltas_classifier.pkl (misma carpeta, entrenado sobre el otro
# benchmark).
mkdir -p EMNgly/checkpoints
curl -L -o EMNgly/checkpoints/N-GlyDE.pickle \
  "https://drive.usercontent.google.com/download?id=1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw&export=download&confirm=t"
```

Kinase Library (corroboracion OPCIONAL e informativa de especificidad de
quinasa para fosforilacion, AMBOS caminos -- FASTA y PDB, analisis de
coherencia biologica 2026-08-07 punto 5, ver seccion "Alcance e
interpretacion" arriba -- NUNCA un motor de consenso,
`Settings.KINASE_LIBRARY_ENABLED` degrada solo con un aviso si no esta
instalado). Entorno DEDICADO: `numpy~=1.26.4`/`pandas~=2.2.3` que fija el
propio paquete son incompatibles con las versiones fijadas de
requirements.txt de este pipeline:

```bash
conda create -n kinase_library python=3.10 -y
conda run -n kinase_library pip install kinase-library
export KINASE_LIBRARY_PYTHON_BIN=$(conda run -n kinase_library which python)
```

> Los pesos del SVM se generaron con `scikit-learn==1.1.1` -- el
> `environment.yml` del repo pinnea 1.5.1, pero cargar los pickles con una
> version distinta a la de entrenamiento dispara `InconsistentVersionWarning`
> (no garantizado). Pinnea `scikit-learn==1.1.1` en el venv, no lo que pide
> el `environment.yml` del repo.

Ver `STATUS.md` para el detalle completo de que ya se verifico (todos los
repos/envs clonados y probados en esta maquina, mas el recurso externo de
StackGlyEmbed, mas los 2 go/no-go checks de EMNGly -- ambos PASARON
2026-08-07) y que falta.

## Licencias

Este repositorio (codigo propio, ver `LICENSE`): CC BY-NC 4.0, uso no
comercial -- consistente con la dependencia real mas restrictiva (ver abajo).

- **DeepMVP**: GPL-3.0 (declarada en `DeepMVP/LICENSE` del repo original -- corregido
  2026-08-03, una version anterior de este README decia MIT por error, nunca
  verificado contra el archivo real hasta ahora).
- **DeepPTMPred**: el repo (`github.com/kuikui-wang/DeepPTMPred`, vendorizado
  en `DeepPTMPred/`) no declara licencia propia, pero el paper asociado es
  CC BY-NC 4.0 (Oxford University Press). Se contacto directamente a los
  autores de correspondencia (Yong Liu, Junwen Wang) para confirmar que el
  codigo sigue los mismos terminos. **Junwen Wang confirmo por email el
  2026-07-29**: *"I confirm that the GitHub code follows the same CC BY-NC
  terms."* Uso no comercial (investigacion/TFG + integracion futura como
  plugin de Scipion por el CNB, institucion publica) encaja sin problema
  dentro de CC BY-NC. Detalle completo en `STATUS.md`.
- **MeToken**: MIT (declarada en `MeToken/LICENSE` del repo original).
- **EMNGly**: sin LICENSE declarado en ningun repo (`StellaHxy/EMNgly` ni su
  duplicado `Xiaoyang878/EMNgly`) -- el paper SI es CC BY 4.0 real
  (confirmado via el XML de PMC, PMC10627407, no solo asumido) pero eso NO
  cubre el codigo. Correo redactado 2026-08-07 a los autores de
  correspondencia reales (Yaojun Wang, `wangyaojun@cau.edu.cn`, China
  Agricultural University; Shiwei Sun, `dwsun@ict.ac.cn`, ICT-CAS -- ambos
  verificados via el mismo XML de PMC, no la pagina de Oxford Academic,
  bloqueada por Cloudflare a fetch directo), envio programado para el
  2026-08-10 (lunes). Mismo patron que funciono con DeepPTMPred -- NO
  bloqueante (a diferencia de CoNglyPred, los pesos ya son publicos y
  descargables sin depender de esa respuesta). MIF (vendorizado dentro de
  EMNgly, `model/MIF/`) es de Microsoft (`microsoft/protein-sequence-models`),
  licencia BSD-2 permisiva verificada en el repo oficial -- la copia de
  EMNgly perdio su LICENSE al vendorizarlo.
- **PyRosetta** (feature de SASA por residuo dentro de DeepPTMPred): licencia
  academica/no-comercial de RosettaCommons, ya cubierta por el uso de
  investigacion/TFG de este proyecto.

## Decisiones de arquitectura

Documentadas en el vault, no en este repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
